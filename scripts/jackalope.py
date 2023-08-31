from modules.Gfa import iterate_gfa_edges,iterate_gfa_nodes
from modules.Gaf import GafElement, iter_gaf_alignments
from modules.IncrementalIdMap import IncrementalIdMap

from collections import defaultdict
import sys

import matplotlib

from PyQt5.QtCore import Qt, QLineF, QEvent
from PyQt5.QtGui import QBrush, QPainter, QPen, QColor
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
    QHBoxLayout,
    QPushButton,
    QSlider,
    QGridLayout,
    QVBoxLayout,
    QWidget,
    QSpacerItem,
    QSizePolicy
)
import networkx


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


class Window(QWidget):
    def __init__(self):
        super().__init__()

        self.line_width = 3

        self.scene_top = QGraphicsScene()
        self.scene_bottom = QGraphicsScene()

        self.view_bottom = QGraphicsView(self.scene_bottom)
        self.view_bottom.setRenderHint(QPainter.Antialiasing)

        self.view_top = QGraphicsView(self.scene_top)
        self.view_top.setRenderHint(QPainter.Antialiasing)

        self.control_panel_left = None
        self.construct_left_control_panel()

        self.control_panel_top = None
        self.construct_top_control_panel()

        # 4x4 grid
        self.grid = QGridLayout(self)

        #     0  1  2  3
        #  0  C  M  M  M
        #  1  C  T  T  T
        #  2  C  T  T  T
        #  3  C  B  B  B
        #  4  C  B  B  B

        # origin is top left AND x and y are swapped in function calls
        self.grid.addLayout(self.control_panel_top, 0, 1, 1, 3)
        self.grid.addLayout(self.control_panel_left, 0, 0, 4, 1)
        self.grid.addWidget(self.view_top,1,1,2,3)
        self.grid.addWidget(self.view_bottom,3,1,2,3)
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

        gfa_path = "/home/ryan/code/jackalope/data/test/simple.gfa"
        self.gfa_path = "/home/ryan/data/test_hapslap/results/bad_tandem/chr20_4265303-4265453/graph_no_empty.gfa"
        # self.gfa_path = None
        self.load_gfa()
        self.draw_graph()

        self.gaf_path = "/home/ryan/data/test_hapslap/results/bad_tandem/chr20_4265303-4265453/reads_vs_graph.gaf"
        # self.gaf_path = None
        self.load_gaf()

        self.bottom_highlight_items = list()

        self.scene_top.selectionChanged.connect(self.on_select_alignment_block)

        self.view_bottom.viewport().installEventFilter(self)

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
        self.scene_top.clear()
        self.alignments = defaultdict(list)
        self.alignment_combobox.blockSignals(True)
        self.alignment_combobox.clear()
        self.alignment_combobox.blockSignals(False)
        self.gaf_query_combobox.blockSignals(True)
        self.gaf_query_combobox.clear()
        self.gaf_query_combobox.blockSignals(False)

    def load_gaf(self):
        self.clear_gaf()

        print("a")
        self.alignments = defaultdict(list)
        for a in iter_gaf_alignments(self.gaf_path):
            self.alignments[a.get_query_name()].append(a)
            print("b")

            for name,reversal in a.get_path():
                if name not in self.qt_nodes:
                    print(a.get_query_name() + " not found in graph")
                    d = OkPopup("ERROR", "GAF contains node names which do not exist in GFA, or GFA is not loaded. \n\n"
                                         "Please open compatible GFA first.")
                    d.exec()
                    self.clear_gaf()

                    return

        print("c")

        self.gaf_query_combobox.addItems(sorted(self.alignments.keys()))

        for query_name in self.alignments.keys():
            self.alignments[query_name] = sorted(self.alignments[query_name], key=lambda x: x.get_query_midpoint())

        print("d")

        # Initialize the menu with whichever query is first
        self.on_select_gaf_query()

    def on_select_gaf_query(self):
        query_name = str(self.gaf_query_combobox.currentText())
        print(query_name)
        print("e")

        self.alignment_combobox.blockSignals(True)
        self.alignment_combobox.clear()

        self.alignment_combobox.addItem("all")
        for i in range(len(self.alignments[query_name])):
            self.alignment_combobox.addItem(str(i))

        self.alignment_combobox.blockSignals(False)

        print("f")

        self.plot_alignments(query_name)
        print("f4")

        self.on_select_alignment()

    def plot_alignments(self, query_name):
        self.scene_top.clear()
        print("f1")

        scale = 1000

        rect = QGraphicsRectItem(0, 0, scale, 20)
        brush = QBrush(Qt.gray)
        rect.setBrush(brush)

        # # Define the pen (line)
        # pen = QPen(Qt.black)
        # pen.setWidth(10)
        # rect.setPen(pen)

        self.scene_top.addItem(rect)
        print("f2")

        for i,alignment in enumerate(self.alignments[query_name]):
            a = alignment.get_query_start()
            b = alignment.get_query_stop()
            l = alignment.get_query_length()

            x1 = float(a)/float(l)*scale
            x2 = float(b)/float(l)*scale
            y = (i+1)*25
            w = x2 - x1

            print(x1, y, w, 20)

            rect2 = QGraphicsRectItem(x1, y, w, 20)
            brush = QBrush(Qt.blue)
            rect2.setBrush(brush)

            print(rect2)

            rect2.instance_item = i
            rect2.setFlag(QGraphicsItem.ItemIsSelectable)

            self.scene_top.addItem(rect2)

        print("f3")

    def color_alignment(self, alignment: GafElement):
        path = alignment.get_path()
        print(path)

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
        print(path)

        for p,[name,reversal] in enumerate(path):
            if name in self.qt_nodes:
                line = self.qt_nodes[name]

                pen = line.pen()
                pen.setWidth(self.line_width+1)
                pen.setColor(Qt.black)

                item = self.scene_bottom.addLine(line.line(), pen)
                item.setZValue(line.zValue()-1)
                self.bottom_highlight_items.append(item)

            else:
                raise Exception("ERROR: bad GAF node")

    def on_select_alignment_block(self):
        for item in self.bottom_highlight_items:
            self.scene_bottom.removeItem(item)

        self.bottom_highlight_items = list()

        items = self.scene_top.selectedItems()

        for item in items:
            if item is not None:
                alignment_index = item.instance_item

                query_name = str(self.gaf_query_combobox.currentText())
                a = self.alignments[query_name][alignment_index]

                self.highlight_alignment(a)

                print(alignment_index)

    def on_select_alignment(self):
        # First reset the colors
        for node in self.qt_nodes.values():
            pen = QPen(Qt.gray)
            pen.setWidth(self.line_width)
            node.setPen(pen)

        query_name = str(self.gaf_query_combobox.currentText())
        selection = self.alignment_combobox.currentText()

        print("g")

        if selection == "all":
            for a,alignment in enumerate(self.alignments[query_name]):
                self.color_alignment(alignment)

        else:
            alignment_index = int(selection)
            alignment = self.alignments[query_name][alignment_index]
            self.color_alignment(alignment)

        print("h")

    def redraw_graph(self):
        if self.gfa_path is not None and len(self.sequences) > 0:
            self.scene_bottom.clear()
            self.draw_graph()
            self.on_select_gaf_query()

    def construct_left_control_panel(self):
        self.control_panel_left = QVBoxLayout()

        up = QPushButton("Open GFA (graph)")
        up.clicked.connect(self.open_gfa)
        self.control_panel_left.addWidget(up)

        up = QPushButton("Open GAF (alignments)")
        up.clicked.connect(self.open_gaf)
        self.control_panel_left.addWidget(up)

        up = QPushButton("Redraw graph")
        up.clicked.connect(self.redraw_graph)
        self.control_panel_left.addWidget(up)

        self.control_panel_left.setAlignment(Qt.AlignTop)
        self.control_panel_left.setSpacing(2)
        self.control_panel_left.addStretch(1)

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

        return super().eventFilter(source,event)

    def draw_graph(self):
        interval_size = 500
        scale = 100

        seed_graph = networkx.Graph()
        graph = networkx.Graph()

        for name,sequence in self.sequences.items():
            seed_graph.add_node(name)

            n = len(sequence) // interval_size + 1
            w = 10

            graph.add_node(name + "_left")
            graph.add_node(name + "_right")
            # graph.add_node(sequence.name + "_top")
            # graph.add_node(sequence.name + "_bottom")

            graph.add_edge(name + "_left",name + "_right", weight=-0.01)
            # graph.add_edge(sequence.name + "_top",sequence.name + "_bottom", weight=-1)
            # graph.add_edge(sequence.name + "_left",sequence.name + "_bottom", weight=-1)
            # graph.add_edge(sequence.name + "_right",sequence.name + "_bottom", weight=-1)
            # graph.add_edge(sequence.name + "_left",sequence.name + "_top", weight=-1)
            # graph.add_edge(sequence.name + "_right",sequence.name + "_top", weight=-1)

            for i in range(n):
                top_name = name + "_" + str(i) + "_top"
                bottom_name = name + "_" + str(i) + "_bottom"

                graph.add_node(top_name)
                graph.add_node(bottom_name)
                graph.add_edge(top_name,bottom_name)

                if i > 0:
                    prev_top_name = name + "_" + str(i-1) + "_top"
                    prev_bottom_name = name + "_" + str(i-1) + "_bottom"

                    graph.add_edge(top_name,prev_bottom_name)
                    graph.add_edge(top_name,prev_top_name)
                    # graph.add_edge(bottom_name,prev_top_name)
                    graph.add_edge(bottom_name,prev_bottom_name)

                if i == 0:
                    graph.add_edge(name + "_left", top_name, weight=w)
                    graph.add_edge(name + "_left", bottom_name, weight=w)

                if i == n-1:
                    graph.add_edge(name + "_right", top_name, weight=w)
                    graph.add_edge(name + "_right", bottom_name, weight=w)

        for edge in self.edges:
            l_a = len(self.sequences[edge.name_a])
            l_b = len(self.sequences[edge.name_b])
            seed_graph.add_edge(edge.name_a, edge.name_b, weight=1/(max(l_a,l_b)))

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

            graph.add_edge(a,b)

        seed_layout = networkx.spring_layout(seed_graph, weight='weight', k=1/len(self.sequences))
        layout = networkx.spring_layout(graph, weight='weight', pos=seed_layout, iterations=2000, k=1/graph.number_of_nodes())

        # for a,b in graph.edges():
        #     pen = QPen(Qt.red)
        #     x_a,y_a = layout[a]
        #     x_b,y_b = layout[b]
        #
        #     ellipse = QGraphicsEllipseItem(x_a*scale,y_a*scale,1,1)
        #     self.scene_bottom.addItem(ellipse)
        #     ellipse = QGraphicsEllipseItem(x_b*scale,y_b*scale,1,1)
        #     self.scene_bottom.addItem(ellipse)
        #
        #     self.scene_bottom.addLine(QLineF(x_a*scale, y_a*scale, x_b*scale, y_b*scale), pen)

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

            item = self.scene_bottom.addLine(QLineF(x_a*scale, y_a*scale, x_b*scale, y_b*scale))

        for name in self.sequences.keys():
            x_a,y_a = layout[name + "_left"]
            x_b,y_b = layout[name + "_right"]
            pen = QPen(Qt.gray)
            pen.setWidth(self.line_width)

            item = self.scene_bottom.addLine(QLineF(x_a*scale, y_a*scale, x_b*scale, y_b*scale), pen)
            item.setFlag(QGraphicsItem.ItemIsSelectable)

            self.qt_nodes[name] = item


def main():
    app = QApplication(sys.argv)
    w = Window()
    w.show()
    app.exec()


if __name__ == "__main__":
    main()
