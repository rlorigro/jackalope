import os.path
import random

from modules.Gfa import iterate_gfa_edges,iterate_gfa_nodes
from modules.Gaf import GafElement, iter_gaf_alignments
from modules.IncrementalIdMap import IncrementalIdMap
from cugraph import Graph,force_atlas2
from cudf import DataFrame
from collections import defaultdict
import numpy
import math
import sys

import matplotlib

from PyQt5.QtCore import Qt, QLineF, QEvent, QUrl
from PyQt5.QtGui import QBrush, QPainter, QPen, QColor, QPainterPath
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import (
    QDialog,
    QLabel,
    QDialogButtonBox,
    QFileDialog,
    QComboBox,
    QApplication,
    QGraphicsEllipseItem,
    QGraphicsItem,
    QGraphicsRectItem,
    QGraphicsScene,
    QGraphicsView,
    QPlainTextEdit,
    QHBoxLayout,
    QPushButton,
    QSlider,
    QGridLayout,
    QVBoxLayout,
    QWidget,
    QSpacerItem,
    QSizePolicy
)


def get_midpoint(a,b):
    x = float(a[0] + b[0])/2.0
    y = float(a[1] + b[1])/2.0

    return x,y


class OkPopup(QDialog):
    def __init__(self, title, message):
        super().__init__()

        self.setWindowTitle(title)

        QBtn = QDialogButtonBox.Ok

        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)

        self.layout = QVBoxLayout()
        message = QLabel(message)
        self.layout.addWidget(message)
        self.layout.addWidget(self.buttonBox)
        self.setLayout(self.layout)


class GafStatsPopup(QDialog):
    def __init__(self, title, alignment: GafElement):
        super().__init__()

        self.setWindowTitle(title)

        QBtn = QDialogButtonBox.Ok

        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)

        self.layout = QVBoxLayout()

        self.cigar_label = QLabel("Cigar")
        self.cigar_box = QPlainTextEdit()

        self.path_label = QLabel("Path")
        self.path_box = QPlainTextEdit()

        cigar_string = ""

        n_match = 0
        n_mismatch = 0
        n_delete = 0
        n_insert = 0
        for c,length in alignment.get_cigar():
            cigar_string += c
            cigar_string += str(length)

            if c == "=":
                n_match += length
            if c == "X":
                n_mismatch += length
            if c == "D":
                n_delete += length
            if c == "I":
                n_insert += length

        identity = float(n_match) / float(n_match + n_mismatch + n_insert + n_delete)

        self.stats_text = QLabel(
            "n_match:\t" + str(n_match) + '\n'
            "n_mismatch:\t" + str(n_mismatch) + '\n'
            "n_delete:\t" + str(n_delete) + '\n'
            "n_insert:\t\t" + str(n_insert) + '\n'
            "identity:\t\t" + "%.5f" % identity + '\n'
        )

        self.path_box.setPlainText(alignment.get_path_string())
        self.cigar_box.setPlainText(cigar_string)

        self.layout.addWidget(self.path_label)
        self.layout.addWidget(self.path_box)

        self.layout.addWidget(self.cigar_label)
        self.layout.addWidget(self.cigar_box)

        self.layout.addWidget(self.stats_text)

        self.layout.addWidget(self.buttonBox)
        self.setLayout(self.layout)


class Window(QWidget):
    def __init__(self):
        super().__init__()

        self.line_width = 2
        self.highlight_width = 1

        self.scene_middle = QGraphicsScene()
        self.scene_bottom = QGraphicsScene()

        # self.project_dir = os.path.dirname(os.path.dirname(__file__))
        # self.igv_html = os.path.join(self.project_dir, "html/igv.html")
        #
        # print(self.igv_html)
        #
        # self.view_top = QWebEngineView()
        # self.view_top.resize(800,300)
        # # self.view_top.load(QUrl("http://www.google.com/"))
        # self.view_top.load(QUrl.fromLocalFile(self.igv_html))
        # self.view_top.show()
        # self.view_top.loadFinished.connect(self.load_finished)

        self.view_bottom = QGraphicsView(self.scene_bottom)
        self.view_bottom.setRenderHint(QPainter.Antialiasing)

        self.view_middle = QGraphicsView(self.scene_middle)
        self.view_middle.setRenderHint(QPainter.Antialiasing)

        self.control_panel_left = None
        self.construct_left_control_panel()

        self.control_panel_top = None
        self.construct_top_control_panel()

        # 4x4 grid
        self.grid = QGridLayout(self)

        #     0  1  2  3
        #  0  C  X  X  X
        #  1  C  T  T  T
        #  2  C  T  T  T
        #  3  C  M  M  M
        #  4  C  B  B  B
        #  5  C  B  B  B

        # origin is top left AND x and y are swapped in function calls
        self.grid.addLayout(self.control_panel_top, 0,1,1,3)
        self.grid.addLayout(self.control_panel_left, 0,0,4,1)
        # self.grid.addWidget(self.view_top, 1,1,2,3)
        self.grid.addWidget(self.view_middle, 1,1,1,3)
        self.grid.addWidget(self.view_bottom, 3,1,2,3)
        self.setLayout(self.grid)

        self.gaf_query_combobox = QComboBox()
        self.gaf_query_combobox.currentIndexChanged.connect(self.on_select_gaf_query)
        self.control_panel_top.addWidget(self.gaf_query_combobox)

        self.alignment_combobox = QComboBox()
        self.alignment_combobox.currentIndexChanged.connect(self.on_select_alignment)
        self.control_panel_top.addWidget(self.alignment_combobox)

        self.colormap = matplotlib.colormaps['jet']

        # Graph data structures
        self.sequences = dict()
        self.edges = list()
        self.name_to_edge = defaultdict(list)
        self.qt_nodes = dict()

        # Alignment data
        self.alignments = defaultdict(list)

        # self.gfa_path = "/home/ryan/data/test_hapslap/results/bad_tandem/chr20_4265303-4265453/graph_no_empty.gfa"
        self.gfa_path = "/home/ryan/data/test_hapslap/evaluation/competitors/all/regional/chr20_65953277-65953889/sniffles_raw/variant_graph.gfa"
        # self.gfa_path = None
        self.load_gfa()
        self.draw_graph()

        # self.gaf_path = "/home/ryan/data/test_hapslap/results/bad_tandem/chr20_4265303-4265453/reads_vs_graph.gaf"
        self.gaf_path = "/home/ryan/data/test_hapslap/evaluation/competitors/all/regional/chr20_65953277-65953889/sniffles_raw/reads_vs_graph.gaf"
        # self.gaf_path = None
        self.load_gaf()

        self.bottom_highlight_items = list()

        self.scene_middle.selectionChanged.connect(self.on_select_alignment_block)

        self.view_bottom.viewport().installEventFilter(self)

    def load_finished(self):
        print("done")

    def open_gfa(self):
        # https://pythonspot.com/pyqt5-file-dialog/
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        filename, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","GFA Files (*.gfa)", options=options)

        if filename != '':
            self.gfa_path = filename
            self.load_gfa()
            self.draw_graph()

    def open_gaf(self):
        # https://pythonspot.com/pyqt5-file-dialog/
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        filename, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","GAF Files (*.gaf)", options=options)

        if filename != '':
            self.gaf_path = filename
            self.load_gaf()

    def load_gfa(self):
        # This will be called repeatedly in some cases, so reset the relevant datastructures
        self.scene_bottom.clear()
        self.bottom_highlight_items = list()

        self.qt_nodes = dict()
        self.sequences = dict()
        self.edges = list()
        self.name_to_edge = defaultdict(list)

        for sequence in iterate_gfa_nodes(self.gfa_path):
            self.sequences[sequence.name] = sequence.sequence

        for edge in iterate_gfa_edges(self.gfa_path):
            self.edges.append(edge)
            self.name_to_edge[edge.name_a].append(len(self.edges) - 1)
            self.name_to_edge[edge.name_b].append(len(self.edges) - 1)

        self.clear_gaf()

    def clear_gaf(self):
        self.scene_middle.clear()
        self.alignments = defaultdict(list)
        self.alignment_combobox.blockSignals(True)
        self.alignment_combobox.clear()
        self.alignment_combobox.blockSignals(False)
        self.gaf_query_combobox.blockSignals(True)
        self.gaf_query_combobox.clear()
        self.gaf_query_combobox.blockSignals(False)

    def load_gaf(self):
        self.clear_gaf()

        self.alignments = defaultdict(list)
        for a in iter_gaf_alignments(self.gaf_path):
            self.alignments[a.get_query_name()].append(a)

            for name,reversal in a.get_path():
                if name not in self.qt_nodes:
                    print(a.get_query_name() + " not found in graph")
                    d = OkPopup("ERROR", "GAF contains node names which do not exist in GFA, or GFA is not loaded. \n\n"
                                         "Please open compatible GFA first.")
                    d.exec()
                    self.clear_gaf()

                    return

        self.gaf_query_combobox.addItems(sorted(self.alignments.keys()))

        for query_name in self.alignments.keys():
            self.alignments[query_name] = sorted(self.alignments[query_name], key=lambda x: x.get_query_midpoint())

        # Initialize the menu with whichever query is first
        self.on_select_gaf_query()

    def on_select_gaf_query(self):
        query_name = str(self.gaf_query_combobox.currentText())

        self.alignment_combobox.blockSignals(True)
        self.alignment_combobox.clear()

        self.alignment_combobox.addItem("all")
        for i in range(len(self.alignments[query_name])):
            self.alignment_combobox.addItem(str(i))

        self.alignment_combobox.blockSignals(False)

        self.plot_alignments(query_name)
        self.on_select_alignment()

    def plot_alignments(self, query_name):
        self.scene_middle.clear()

        scale = 1000

        rect = QGraphicsRectItem(0, 0, scale, 20)
        brush = QBrush(Qt.gray)
        rect.setBrush(brush)

        # # Define the pen (line)
        # pen = QPen(Qt.black)
        # pen.setWidth(10)
        # rect.setPen(pen)

        self.scene_middle.addItem(rect)

        for i,alignment in enumerate(self.alignments[query_name]):
            a = alignment.get_query_start()
            b = alignment.get_query_stop()
            l = alignment.get_query_length()

            x1 = float(a)/float(l)*scale
            x2 = float(b)/float(l)*scale
            y = (i+1)*25
            w = x2 - x1

            rect2 = QGraphicsRectItem(x1, y, w, 20)
            brush = QBrush(Qt.blue)
            rect2.setBrush(brush)

            rect2.instance_item = i
            rect2.setFlag(QGraphicsItem.ItemIsSelectable)

            self.scene_middle.addItem(rect2)

    def color_alignment(self, alignment: GafElement):
        path = alignment.get_path()

        for p,[name,reversal] in enumerate(path):
            color = self.colormap(float(p)/float(len(path)))

            color = [int(round(255*c)) for c in color]

            pen = QPen(QColor.fromRgb(color[0], color[1], color[2]))
            pen.setWidth(self.line_width)

            if name in self.qt_nodes:
                self.qt_nodes[name].setPen(pen)
            else:
                raise Exception("ERROR: bad GAF node")

    def highlight_alignment(self, alignment: GafElement):
        path = alignment.get_path()

        for p,[name,reversal] in enumerate(path):
            if name in self.qt_nodes:
                path = self.qt_nodes[name]

                pen = path.pen()
                pen.setWidth(self.line_width+self.highlight_width)
                pen.setColor(Qt.black)

                item = self.scene_bottom.addPath(path.path(), pen)
                item.setZValue(path.zValue()-1)
                self.bottom_highlight_items.append(item)

            else:
                raise Exception("ERROR: bad GAF node")

    def on_select_alignment_block(self):
        for item in self.bottom_highlight_items:
            self.scene_bottom.removeItem(item)

        self.bottom_highlight_items = list()

        items = self.scene_middle.selectedItems()

        for item in items:
            if item is not None:
                alignment_index = item.instance_item

                query_name = str(self.gaf_query_combobox.currentText())
                a = self.alignments[query_name][alignment_index]

                self.highlight_alignment(a)

    def on_select_alignment(self):
        # First reset the colors
        for node in self.qt_nodes.values():
            pen = QPen(Qt.gray)
            pen.setWidth(self.line_width)
            node.setPen(pen)

        query_name = str(self.gaf_query_combobox.currentText())
        selection = self.alignment_combobox.currentText()

        if selection == "all":
            for a,alignment in enumerate(self.alignments[query_name]):
                self.color_alignment(alignment)

        else:
            alignment_index = int(selection)
            alignment = self.alignments[query_name][alignment_index]
            self.color_alignment(alignment)

    def redraw_graph(self):
        if self.gfa_path is not None and len(self.sequences) > 0:
            self.scene_bottom.clear()
            self.draw_graph()
            self.on_select_gaf_query()

    def construct_left_control_panel(self):
        self.control_panel_left = QVBoxLayout()

        button = QPushButton("Open GFA (graph)")
        button.clicked.connect(self.open_gfa)
        self.control_panel_left.addWidget(button)

        button = QPushButton("Open GAF (alignments)")
        button.clicked.connect(self.open_gaf)
        self.control_panel_left.addWidget(button)

        button = QPushButton("Redraw graph")
        button.clicked.connect(self.redraw_graph)
        self.control_panel_left.addWidget(button)

        button = QPushButton("Show alignment details")
        button.clicked.connect(self.show_alignment_details)
        self.control_panel_left.addWidget(button)

        self.control_panel_left.setAlignment(Qt.AlignTop)
        self.control_panel_left.setSpacing(2)
        self.control_panel_left.addStretch(1)

    def show_alignment_details(self):
        items = self.scene_middle.selectedItems()

        if len(items) == 0:
            d = OkPopup("ERROR", "No alignment blocks selected, please select an alignment first")
            d.exec()

        elif len(items) == 1:
            alignment_index = items[0].instance_item

            query_name = str(self.gaf_query_combobox.currentText())
            a = self.alignments[query_name][alignment_index]

            d = GafStatsPopup("Alignment details", alignment=a)
            d.exec()

        else:
            d = OkPopup("ERROR", "Too many alignments selected, please select one at a time")
            d.exec()

    def construct_top_control_panel(self):
        self.control_panel_top = QVBoxLayout()

        self.control_panel_top.setAlignment(Qt.AlignTop)
        self.control_panel_top.setSpacing(2)
        # self.control_panel_top.addStretch(1)

    def eventFilter(self, source, event):
        # Combined from 2 sources:
        # https://stackoverflow.com/questions/63634323/zooming-image-on-qgraphicsscene-with-mouse-wheel-in-pyqt5
        # https://stackoverflow.com/questions/58965209/zoom-on-mouse-position-qgraphicsview

        if (source == self.view_bottom.viewport() and event.type() == QEvent.Wheel) and event.modifiers() == Qt.ControlModifier:
            factor = 1.1
            if event.angleDelta().y() < 0:
                factor = 0.9

            view_pos = event.pos()
            scene_pos = self.view_bottom.mapToScene(view_pos)
            self.view_bottom.centerOn(scene_pos)
            self.view_bottom.scale(factor, factor)
            delta = self.view_bottom.mapToScene(view_pos) - self.view_bottom.mapToScene(self.view_bottom.viewport().rect().center())
            self.view_bottom.centerOn(scene_pos - delta)            # do not propagate the event to the scroll area scrollbars
            return True

        elif event.type() == QEvent.GraphicsSceneMousePress:
            pass

        elif source == self.view_middle.viewport() and event.type() == event.Type.ContextMenu:
            self.view_middle.contextMenuEvent(event)

            return event.isAccepted()

        return super().eventFilter(source,event)

    def build_path(self, points):
        path = QPainterPath()

        for i,[x,y] in enumerate(points):
            if i == 0:
                path.moveTo(x,y)
            else:
                path.lineTo(x,y)

        return path

    def draw_graph(self):
        # interval_size = 500
        scale = 0.03

        seed_graph = Graph()
        seed_edges_as_ids = list()
        seed_id_map = IncrementalIdMap()

        graph = Graph()
        edges_as_ids = list()
        id_map = IncrementalIdMap()

        total_length = 0

        for sequence in self.sequences.values():
            total_length += len(sequence)

        interval_size = total_length//100 + 1

        subnodes = defaultdict(list)

        for name,sequence in self.sequences.items():
            n = max(3,len(sequence) // interval_size)
            w = 1

            name_seed = name
            seed_id_map.add(name)

            name_left = name + "_left"
            name_right = name + "_right"

            id_left = id_map.add(name_left)
            id_right = id_map.add(name_right)

            # edges_as_ids.append((id_left,id_right,-0.01))

            for i in range(n):
                name_top = name + "_" + str(i) + "_top"
                name_bottom = name + "_" + str(i) + "_bottom"

                id_top = id_map.add(name_top)
                id_bottom = id_map.add(name_bottom)

                # Keep track of which middle pairs of nodes constitute the scaffolding for a given node
                subnodes[name].append([name_top,name_bottom])

                edges_as_ids.append((id_top,id_bottom,w))

                if i > 0:
                    prev_id_top = id_map.get_id(name + "_" + str(i-1) + "_top")
                    prev_id_bottom = id_map.get_id(name + "_" + str(i-1) + "_bottom")

                    edges_as_ids.append((id_top, prev_id_top, w))
                    edges_as_ids.append((id_bottom, prev_id_bottom, w))
                    edges_as_ids.append((id_top, prev_id_bottom, w))
                    edges_as_ids.append((id_bottom, prev_id_top, w))

                if i == 0:
                    edges_as_ids.append((id_left, id_top, w))
                    edges_as_ids.append((id_left, id_bottom, w))

                if i == n-1:
                    edges_as_ids.append((id_right, id_top, w))
                    edges_as_ids.append((id_right, id_bottom, w))

        for edge in self.edges:
            l_a = len(self.sequences[edge.name_a])
            l_b = len(self.sequences[edge.name_b])

            a = None
            b = None

            if edge.reversal_a:
                a = edge.name_a + "_left"
            else:
                a = edge.name_a + "_right"

            if edge.reversal_b:
                b = edge.name_b + "_right"
            else:
                b = edge.name_b + "_left"

            a_seed = seed_id_map.get_id(edge.name_a)
            b_seed = seed_id_map.get_id(edge.name_b)
            # w_seed = 1 / math.log10(max(len(self.sequences[edge.name_a]), len(self.sequences[edge.name_b])) + 10)

            seed_edges_as_ids.append((a_seed, b_seed, 0.5))

            edges_as_ids.append((id_map.get_id(a), id_map.get_id(b), 0.5))

        seed_edge_df = DataFrame(seed_edges_as_ids)

        edge_df = DataFrame(edges_as_ids)

        print("Loading coarse graph as CuDF...")

        seed_positions = list()
        for id,name in seed_id_map:
            x = random.randint(-50,50)
            y = random.randint(-50,50)
            seed_positions.append((id,x,y))

        seed_df = DataFrame(seed_positions, columns=['vertex','x','y'])

        print()
        print("Coarse graph layout...")

        seed_graph.from_cudf_edgelist(seed_edge_df, source=0, destination=1, weight=2)
        seed_result = force_atlas2(seed_graph, max_iter=1000, pos_list=seed_df, scaling_ratio=8.0, verbose=True)

        print(type(seed_result))
        print(seed_result.shape)

        initial_positions = list()
        prev_name = None
        item = None

        for id,name in id_map:
            seed_name = name.split('_')[0]
            seed_id = seed_id_map.get_id(seed_name)

            if seed_name != prev_name:
                item = seed_result.iloc[seed_id]

            x = item['x']
            y = item['y']
            # print(seed_name, item)

            # print(seed_name, x, y)

            initial_positions.append((id,x,y))

            prev_name = seed_name

        print("Loading full graph as CuDF...")
        initial_positions_df = DataFrame(initial_positions, columns=['vertex','x','y'])

        print()
        print("Graph layout...")

        print(edge_df.shape)
        # print(df[0])
        # exit()
        # adjacency_matrix = numpy.zeros([len(id_map)]*2, dtype=float)

        # for a,b,w in edges_as_ids:
        #     adjacency_matrix[a,b] = w
        #     adjacency_matrix[b,a] = w

        graph.from_cudf_edgelist(edge_df, source=0, destination=1, weight=2)

        result = force_atlas2(graph, pos_list=initial_positions_df, max_iter=4000, jitter_tolerance=0.3, scaling_ratio=4.0, verbose=True, barnes_hut_theta=0.90)

        layout = dict()
        for i in range(len(result["vertex"])):
            id = int(result["vertex"][i])
            name = id_map.get_name(id)

            layout[name] = (scale*result['x'][i], scale*result['y'][i])

        # layout = networkx.spring_layout(graph, weight='weight', pos=seed_layout, iterations=2000, k=1/graph.number_of_nodes())

        # for id_a,id_b,w in edges_as_ids:
        #     a = id_map.get_name(id_a)
        #     b = id_map.get_name(id_b)
        #
        #     pen = QPen(Qt.red)
        #     x_a,y_a = layout[a]
        #     x_b,y_b = layout[b]
        #
        #     ellipse = QGraphicsEllipseItem(x_a,y_a,1,1)
        #     self.scene_bottom.addItem(ellipse)
        #     ellipse = QGraphicsEllipseItem(x_b,y_b,1,1)
        #     self.scene_bottom.addItem(ellipse)
        #
        #     self.scene_bottom.addLine(QLineF(x_a, y_a, x_b, y_b), pen)

        for edge in self.edges:
            x_a = None
            y_a = None
            x_b = None
            y_b = None

            if edge.reversal_a:
                x_a,y_a = layout[edge.name_a + "_left"]
            else:
                x_a,y_a = layout[edge.name_a + "_right"]

            if edge.reversal_b:
                x_b,y_b = layout[edge.name_b + "_right"]
            else:
                x_b,y_b = layout[edge.name_b + "_left"]

            color = QColor(Qt.black)
            color.setAlphaF(0.7)

            pen = QPen(color)
            pen.setWidth(max(1,self.line_width//3))

            item = self.scene_bottom.addLine(QLineF(x_a, y_a, x_b, y_b), pen=pen)

        for name,node_pairs in subnodes.items():
            start = name + "_left"
            stop = name + "_right"

            points = [layout[start]] + [get_midpoint(layout[a],layout[b]) for a,b in node_pairs] + [layout[stop]]

            path = self.build_path(points)

            pen = QPen(Qt.gray | Qt.FlatCap | Qt.BevelJoin)
            pen.setWidth(self.line_width)
            item = self.scene_bottom.addPath(path, pen)

            item.setFlag(QGraphicsItem.ItemIsSelectable)
            self.qt_nodes[name] = item



        # for name in self.sequences.keys():
        #     x_a,y_a = layout[name + "_left"]
        #     x_b,y_b = layout[name + "_right"]
        #     pen = QPen(Qt.gray)
        #     pen.setWidth(self.line_width)
        #
        #     item = self.scene_bottom.addLine(QLineF(x_a*scale, y_a*scale, x_b*scale, y_b*scale), pen)
        #     item.setFlag(QGraphicsItem.ItemIsSelectable)
        #
        #     self.qt_nodes[name] = item


def main():
    sys.argv.append("--disable-web-security")

    app = QApplication(sys.argv)
    w = Window()
    w.resize(1000,600)

    w.show()
    app.exec()


if __name__ == "__main__":
    main()
