import multiprocessing
from collections import defaultdict
import tempfile
import os.path
import random

from modules.Gfa import iterate_gfa_edges,iterate_gfa_nodes
from modules.Gaf import GafElement, iter_gaf_alignments
from modules.IncrementalIdMap import IncrementalIdMap
from modules.Align import run_minigraph,run_panaligner,run_graphaligner

import matplotlib
import networkx

import importlib.util
import numpy
import math
import sys

USE_CUDA = False
if USE_CUDA and importlib.util.find_spec("cugraph") is not None and importlib.util.find_spec("cudf") is not None:
    sys.stderr.write("Found: cugraph/cudf libs\n")
    from cugraph import Graph as cuGraph
    from cugraph import force_atlas2
    from cudf import DataFrame
else:
    sys.stderr.write("NOT found: cugraph/cudf libs\n")


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
    QGraphicsTextItem,
    QGraphicsScene,
    QGraphicsView,
    QLineEdit,
    QCheckBox,
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


def parse_string_as_numeric_positive_integer(s):
    i = None

    if s == "":
        pass
    elif not s.isnumeric():
        d = OkPopup("ERROR", "Must enter numeric integer value")
        d.exec()
    else:
        i = int(s)

        if i < 0:
            d = OkPopup("ERROR", "Must enter POSITIVE numeric integer value")
            d.exec()
            i = None

    return i


def get_midpoint(a,b):
    x = float(a[0] + b[0])/2.0
    y = float(a[1] + b[1])/2.0

    return x,y


class RunAlignmentPopup(QDialog):
    def __init__(self, gfa_path):
        super().__init__()

        self.gfa_path = gfa_path
        self.gaf_path = None

        self.layout = QVBoxLayout()

        self.setWindowTitle("Run new alignment (GraphAligner)")

        args = [
            "-x", "vg"
        ]

        # Node width field
        field_label = QLabel("Arguments:")
        self.args_field = QPlainTextEdit(" ".join(args))

        self.layout.addWidget(field_label)
        self.layout.addWidget(self.args_field)

        # THREADS FIELD
        threads_layout = QHBoxLayout()
        threads_label = QLabel("Threads:")

        self.threads_field = QLineEdit(str(multiprocessing.cpu_count()))
        self.threads_field.textChanged.connect(self.parse_threads)

        # Construct composite row of elements
        threads_layout.addWidget(threads_label)
        threads_layout.addWidget(self.threads_field)

        # Add to dialog box
        self.layout.addLayout(threads_layout)

        temp_dir = tempfile.mkdtemp()

        # INPUT FIELD
        # from unknown location
        input_layout = QHBoxLayout()
        input_label = QLabel("Input Fasta file:")

        self.input_field = QLineEdit("")

        button = QPushButton("Browse")
        button.clicked.connect(self.open_fasta)

        # Construct composite row of elements
        input_layout.addWidget(input_label)
        input_layout.addWidget(self.input_field)
        input_layout.addWidget(button)

        # Add to dialog box
        self.layout.addLayout(input_layout)

        # OUTPUT FIELD
        # to an auto-generated tmp file by default, but allow other options
        output_layout = QHBoxLayout()
        output_label = QLabel("Output location (must exist):")

        self.output_field = QLineEdit(temp_dir)

        button = QPushButton("Browse")
        button.clicked.connect(self.open_directory)

        # Construct composite row of elements
        output_layout.addWidget(output_label)
        output_layout.addWidget(self.output_field)
        output_layout.addWidget(button)

        # Add to dialog box
        self.layout.addLayout(output_layout)

        # Should the alignments replace the current GAF alignments or append it?
        self.replacement_checkbox = QCheckBox("Replace alignments (default=append)")
        self.layout.addWidget(self.replacement_checkbox)

        run_button = QPushButton("Run")
        run_button.clicked.connect(self.run_alignment)
        self.layout.addWidget(run_button)

        QBtn = QDialogButtonBox.Ok

        self.buttonBox = QDialogButtonBox(QBtn)
        self.buttonBox.accepted.connect(self.accept)

        self.layout.addWidget(self.buttonBox)
        self.setLayout(self.layout)

    def parse_threads(self):
        i = parse_string_as_numeric_positive_integer(self.threads_field.text())

        return i

    def open_directory(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontResolveSymlinks
        options |= QFileDialog.ShowDirsOnly
        directory = QFileDialog.getExistingDirectory(self, options=options)

        if directory != '':
            self.output_field.setText(directory)

    def open_fasta(self):
        d = QFileDialog()
        d.setFileMode(QFileDialog.ExistingFile)
        filename, _ = d.getOpenFileName(self)

        if filename != '' and str(filename).endswith(".fasta"):
            self.input_field.setText(filename)

    def run_alignment(self):
        args = self.args_field.toPlainText().strip().split()
        n_threads = self.parse_threads()
        fasta_path = self.input_field.text()
        output_dir = self.output_field.text()

        self.gaf_path = run_graphaligner(
            output_directory=output_dir,
            gfa_path=self.gfa_path,
            fasta_path=fasta_path,
            n_threads=n_threads,
            args_override=args
        )


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
        cigar = None

        try:
            cigar = alignment.get_cigar()
        except Exception as e:
            d = OkPopup("ERROR", str(e))
            d.exec()
            cigar = None

        n_match = 0
        n_mismatch = 0
        n_delete = 0
        n_insert = 0
        identity = 0

        if cigar is not None:
            for c,length in cigar:
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
    def __init__(self, gfa_path=None, gaf_path=None):
        super().__init__()
        self.use_cugraph = False

        self.line_width = 2
        self.highlight_width = 1
        self.length_scale_factor = 100
        self.layout_iterations = 100
        self.min_node_length = 3

        self.scene_middle = QGraphicsScene()
        self.scene_bottom = QGraphicsScene()

        self.view_bottom = QGraphicsView(self.scene_bottom)
        self.view_bottom.resize(600,300)
        self.view_bottom.setRenderHint(QPainter.Antialiasing)

        self.view_middle = QGraphicsView(self.scene_middle)
        self.view_middle.resize(600,300)
        self.view_middle.setRenderHint(QPainter.Antialiasing)

        self.line_width_field = None
        self.length_scale_factor_field = None
        self.control_panel_left = None
        self.construct_left_control_panel()

        self.top_label = None
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
        self.grid.addLayout(self.control_panel_left, 0,0,4,1)
        self.grid.addLayout(self.control_panel_top,  0,1,1,3)
        self.grid.addWidget(self.view_middle,        1,1,1,3)
        self.grid.addWidget(self.view_bottom,        2,1,2,3)
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

        self.gfa_path = gfa_path
        self.gaf_path = gaf_path

        if gfa_path is not None:
            self.load_gfa()
            self.draw_graph()

        if gaf_path is not None:
            self.load_gaf()

        self.bottom_highlight_items = list()

        self.scene_middle.selectionChanged.connect(self.on_select_alignment_block)

        self.view_bottom.viewport().installEventFilter(self)

        for i in range(self.grid.columnCount()):
            w = 0 if i==0 else 8
            self.grid.setColumnMinimumWidth(i,200)
            self.grid.setColumnStretch(i,w)

        for i in range(self.grid.rowCount()):
            w = 0 if i==0 else 8
            h = 50 if i==0 else 200
            self.grid.setRowMinimumHeight(i,h)
            self.grid.setRowStretch(i,w)

    def load_finished(self):
        print("done")

    def open_gfa(self):
        self.clear_highlights()
        self.scene_bottom.clear()

        # https://pythonspot.com/pyqt5-file-dialog/
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        filename, _ = QFileDialog.getOpenFileName(self,"QFileDialog.getOpenFileName()", "","GFA Files (*.gfa)", options=options)

        if filename != '':
            self.gfa_path = filename
            self.load_gfa()
            self.draw_graph()
            self.on_select_gaf_query()

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
        self.clear_highlights()
        self.clear_gaf()

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

    def clear_gaf(self):
        self.scene_middle.clear()
        self.clear_highlights()
        self.alignments = defaultdict(list)
        self.alignment_combobox.blockSignals(True)
        self.alignment_combobox.clear()
        self.alignment_combobox.blockSignals(False)
        self.gaf_query_combobox.blockSignals(True)
        self.gaf_query_combobox.clear()
        self.gaf_query_combobox.blockSignals(False)

    def load_gaf(self, replace=True):
        if replace:
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

                    # TODO: make clearing optional
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

        text_items = list()

        text_offset = 0
        text_target_limit = -10
        scale = 1000

        rect = QGraphicsRectItem(0, 0, scale, 20)
        brush = QBrush(Qt.gray)
        rect.setBrush(brush)

        label = QGraphicsTextItem(query_name)
        label.setPos(text_offset + 0,0)
        text_items.append(label)

        self.scene_middle.addItem(label)

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
            label = QGraphicsTextItem("Alignment " + str(i))
            label.setPos(text_offset,y)
            text_items.append(label)

            brush = QBrush(Qt.blue)
            rect2.setBrush(brush)

            rect2.instance_item = i
            rect2.setFlag(QGraphicsItem.ItemIsSelectable)

            self.scene_middle.addItem(label)
            self.scene_middle.addItem(rect2)

        for item in text_items:
            w = item.boundingRect().width()
            y = item.pos().y()
            item.setPos(text_target_limit - w, y)

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

    def clear_highlights(self):
        for item in self.bottom_highlight_items:
            self.scene_bottom.removeItem(item)
            self.bottom_highlight_items = list()

    def on_select_alignment_block(self):
        self.clear_highlights()

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
            self.clear_highlights()
            self.scene_bottom.clear()
            self.draw_graph()
            self.on_select_gaf_query()

    def run_new_alignment(self):
        if self.gfa_path is None:
            d = OkPopup("ERROR", "Must open a GFA to align to")
            d.exec()

            return

        d = RunAlignmentPopup(self.gfa_path)
        if d.exec_() == QDialog.Accepted:
            print(d.gaf_path)
            self.gaf_path = d.gaf_path

            if d.replacement_checkbox.isChecked():
                self.clear_gaf()
                self.load_gaf()
            else:
                self.load_gaf(replace=False)

    def construct_left_control_panel(self):
        self.control_panel_left = QVBoxLayout()

        button = QPushButton("Open GFA (graph)")
        button.clicked.connect(self.open_gfa)
        self.control_panel_left.addWidget(button)

        button = QPushButton("Open GAF (alignments)")
        button.clicked.connect(self.open_gaf)
        self.control_panel_left.addWidget(button)

        button = QPushButton("Run new alignment")
        button.clicked.connect(self.run_new_alignment)
        self.control_panel_left.addWidget(button)

        button = QPushButton("Redraw graph")
        button.clicked.connect(self.redraw_graph)
        self.control_panel_left.addWidget(button)

        button = QPushButton("Show alignment details")
        button.clicked.connect(self.show_alignment_details)
        self.control_panel_left.addWidget(button)

        # Node width field
        field_layout = QHBoxLayout()
        field_label = QLabel("Node width:")
        self.line_width_field = QLineEdit(str(self.line_width))
        self.line_width_field.textChanged.connect(self.adjust_node_width)
        field_layout.addWidget(field_label)
        field_layout.addWidget(self.line_width_field)
        self.control_panel_left.addLayout(field_layout)

        # Length scale field
        field_layout = QHBoxLayout()
        field_label = QLabel("Target total points:")
        self.length_scale_factor_field = QLineEdit(str(self.length_scale_factor))
        self.length_scale_factor_field.textChanged.connect(self.adjust_length_scale_factor)
        field_layout.addWidget(field_label)
        field_layout.addWidget(self.length_scale_factor_field)
        self.control_panel_left.addLayout(field_layout)

        # Min length field
        field_layout = QHBoxLayout()
        field_label = QLabel("Minimum node length:")
        self.min_node_length_field = QLineEdit(str(self.min_node_length))
        self.min_node_length_field.textChanged.connect(self.adjust_min_node_length)
        field_layout.addWidget(field_label)
        field_layout.addWidget(self.min_node_length_field)
        self.control_panel_left.addLayout(field_layout)

        # field_layout = QHBoxLayout()
        # field_label = QLabel("Graph layout iterations:")
        # self.layout_iterations_field = QLineEdit(str(self.layout_iterations))
        # self.layout_iterations_field.textChanged.connect(self.adjust_layout_iterations)
        # field_layout.addWidget(field_label)
        # field_layout.addWidget(self.layout_iterations_field)
        # self.control_panel_left.addLayout(field_layout)

        self.control_panel_left.setAlignment(Qt.AlignTop)
        self.control_panel_left.setSpacing(2)
        self.control_panel_left.addStretch(1)

    def adjust_node_width(self):
        s = self.line_width_field.text()
        i = self.parse_string_as_numeric_positive_integer(s)

        if i is None:
            return

        self.line_width = i

        for item in self.qt_nodes.values():
            pen = item.pen()
            pen.setWidth(self.line_width)
            item.setPen(pen)

    def adjust_length_scale_factor(self):
        s = self.length_scale_factor_field.text()
        i = self.parse_string_as_numeric_positive_integer(s)

        if i is None:
            return

        self.length_scale_factor = i

    def adjust_min_node_length(self):
        s = self.min_node_length_field.text()
        i = self.parse_string_as_numeric_positive_integer(s)

        if i is None:
            return

        self.min_node_length = i

    def adjust_layout_iterations(self):
        s = self.length_scale_factor_field.text()
        i = self.parse_string_as_numeric_positive_integer(s)

        if i is None:
            return

        self.layout_iterations = i

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
        self.top_label = QLabel("Alignments")
        self.control_panel_top = QVBoxLayout()
        self.control_panel_top.addWidget(self.top_label)
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

    def layout_with_cugraph(self, seed_positions, seed_edges_as_ids, edges_as_ids, id_map, seed_id_map, scale):
        seed_df = DataFrame(seed_positions, columns=['vertex','x','y'])

        seed_edge_df = DataFrame(seed_edges_as_ids)
        edge_df = DataFrame(edges_as_ids)

        seed_graph = cuGraph()
        graph = cuGraph()

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

        graph.from_cudf_edgelist(edge_df, source=0, destination=1, weight=2)

        result = force_atlas2(graph, pos_list=initial_positions_df, max_iter=4000, jitter_tolerance=0.3, scaling_ratio=4.0, verbose=True, barnes_hut_theta=0.90)

        layout = dict()
        for i in range(len(result["vertex"])):
            id = int(result["vertex"][i])
            name = id_map.get_name(id)

            layout[name] = (scale*result['x'][i], scale*result['y'][i])

        return layout

    def layout_with_graphviz(self, edges_as_ids, id_map, scale):
        graph = networkx.Graph()

        print()
        print("Graph layout...")

        graph.add_nodes_from(range(len(id_map)))

        for a,b,w in edges_as_ids:
            graph.add_edge(a,b,weight=w)

        s = random.randint(0,2**16-1)
        result = networkx.drawing.nx_agraph.graphviz_layout(graph, prog="sfdp", args="-Grepulsiveforce=1.2 -Gstart=%d" % s)

        layout = dict()
        for id,item in result.items():
            name = id_map.get_name(id)

            layout[name] = item

        return layout

    def draw_graph(self):
        # interval_size = 500

        seed_edges_as_ids = list()
        seed_id_map = IncrementalIdMap()

        edges_as_ids = list()
        id_map = IncrementalIdMap()

        total_length = 0
        for sequence in self.sequences.values():
            total_length += len(sequence)

        interval_size = (float(total_length)/float(self.length_scale_factor))

        print("total_length:", total_length)
        print("length_scale_factor:", self.length_scale_factor)
        print("interval_size:", interval_size)

        subnodes = defaultdict(list)

        for name,sequence in self.sequences.items():
            n = max(self.min_node_length,int(round(float(len(sequence)) / float(interval_size))))
            w = 1

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

        print("Loading coarse graph as CuDF...")

        seed_positions = list()
        for id,name in seed_id_map:
            x = random.randint(-50,50)
            y = random.randint(-50,50)
            seed_positions.append((id,x,y))

        if self.use_cugraph:
            scale = 0.03
            layout = self.layout_with_cugraph(seed_positions, seed_edges_as_ids, edges_as_ids, id_map, seed_id_map, scale)
        else:
            scale = 10
            layout = self.layout_with_graphviz(edges_as_ids, id_map, scale)

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

        self.scene_bottom.update()


def main():
    sys.argv.append("--disable-web-security")

    app = QApplication(sys.argv)
    w = Window()
    w.resize(1000,600)

    w.show()
    app.exec()


if __name__ == "__main__":
    main()
