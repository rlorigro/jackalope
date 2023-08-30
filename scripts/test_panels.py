from modules.Gfa import iterate_gfa_edges,iterate_gfa_nodes
from collections import defaultdict
from modules.IncrementalIdMap import IncrementalIdMap
import sys

from PyQt5.QtCore import Qt, QLineF, QEvent
from PyQt5.QtGui import QBrush, QPainter, QPen
from PyQt5.QtWidgets import (
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


#
# Adapted from https://www.pythonguis.com/tutorials/pyqt-qgraphics-vector-graphics/
#


class Window(QWidget):
    def __init__(self):
        # gfa_path = "/home/ryan/code/jackalope/data/test/simple.gfa"
        gfa_path = "/home/ryan/data/test_hapslap/results/bad_tandem/chr20_4265303-4265453/graph_no_empty.gfa"

        super().__init__()

        # If we don't set this on creation, we can set it later with .setSceneRect
        self.scene_top = QGraphicsScene()
        self.scene_bottom = QGraphicsScene()

        control_panel = QVBoxLayout()

        up = QPushButton("Up")
        up.clicked.connect(self.up)
        control_panel.addWidget(up)

        down = QPushButton("Down")
        down.clicked.connect(self.down)
        control_panel.addWidget(down)

        control_panel.setAlignment(Qt.AlignTop)
        control_panel.setSpacing(2)
        control_panel.addStretch(1)
        # control_panel.setContentsMargins(0, 0, 0, 0)

        self.view_bottom = QGraphicsView(self.scene_bottom)
        self.view_bottom.setRenderHint(QPainter.Antialiasing)

        self.view_top = QGraphicsView(self.scene_top)
        self.view_top.setRenderHint(QPainter.Antialiasing)

        # 4x4 grid
        grid = QGridLayout(self)

        #     0  1  2  3
        #  0  C  T  T  T
        #  1  C  T  T  T
        #  2  C  B  B  B
        #  3  C  B  B  B
        grid.addLayout(control_panel,0,0,4,1)
        grid.addWidget(self.view_top,0,1,2,3)
        grid.addWidget(self.view_bottom,2,1,2,3)

        self.setLayout(grid)

        self.draw_graph(gfa_path)

        self.view_bottom.viewport().installEventFilter(self)

    def eventFilter(self, source, event):
        if (source == self.view_bottom.viewport() and event.type() == QEvent.Wheel) and event.modifiers() == Qt.ControlModifier:
            if event.angleDelta().y() > 0:
                scale = 1.25
            else:
                scale = .8

            self.view_bottom.scale(scale, scale)
            # do not propagate the event to the scroll area scrollbars
            return True

        elif event.type() == QEvent.GraphicsSceneMousePress:
            pass

        return super().eventFilter(source,event)

    def draw_graph(self, gfa_path):
        interval_size = 500
        scale = 100

        seed_graph = networkx.Graph()
        graph = networkx.Graph()

        lengths = dict()
        nodes = list()

        for sequence in iterate_gfa_nodes(gfa_path):
            seed_graph.add_node(sequence.name)

            nodes.append(sequence.name)
            lengths[sequence.name] = len(sequence)

            n = len(sequence) // interval_size + 1
            w = 10

            graph.add_node(sequence.name + "_left")
            graph.add_node(sequence.name + "_right")
            # graph.add_node(sequence.name + "_top")
            # graph.add_node(sequence.name + "_bottom")

            graph.add_edge(sequence.name + "_left",sequence.name + "_right", weight=-0.01)
            # graph.add_edge(sequence.name + "_top",sequence.name + "_bottom", weight=-1)
            # graph.add_edge(sequence.name + "_left",sequence.name + "_bottom", weight=-1)
            # graph.add_edge(sequence.name + "_right",sequence.name + "_bottom", weight=-1)
            # graph.add_edge(sequence.name + "_left",sequence.name + "_top", weight=-1)
            # graph.add_edge(sequence.name + "_right",sequence.name + "_top", weight=-1)

            for i in range(n):
                top_name = sequence.name + "_" + str(i) + "_top"
                bottom_name = sequence.name + "_" + str(i) + "_bottom"

                graph.add_node(top_name)
                graph.add_node(bottom_name)
                graph.add_edge(top_name,bottom_name)

                if i > 0:
                    prev_top_name = sequence.name + "_" + str(i-1) + "_top"
                    prev_bottom_name = sequence.name + "_" + str(i-1) + "_bottom"

                    graph.add_edge(top_name,prev_bottom_name)
                    graph.add_edge(top_name,prev_top_name)
                    # graph.add_edge(bottom_name,prev_top_name)
                    graph.add_edge(bottom_name,prev_bottom_name)

                if i == 0:
                    graph.add_edge(sequence.name + "_left", top_name, weight=w)
                    graph.add_edge(sequence.name + "_left", bottom_name, weight=w)

                if i == n-1:
                    graph.add_edge(sequence.name + "_right", top_name, weight=w)
                    graph.add_edge(sequence.name + "_right", bottom_name, weight=w)

        edges = list()
        for edge in iterate_gfa_edges(gfa_path):
            seed_graph.add_edge(edge.name_a, edge.name_b, weight=1/(max(lengths[edge.name_a], lengths[edge.name_b])))

            print(edge)
            edges.append(edge)

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

            print("adding edge:", a, b)
            graph.add_edge(a,b)

        seed_layout = networkx.spring_layout(seed_graph, weight='weight')
        layout = networkx.spring_layout(graph, weight='weight', pos=seed_layout, iterations=2000, k=1/graph.number_of_nodes())

        for item in self.scene_bottom.items():
            item.setFlag(QGraphicsItem.ItemIsMovable)
            item.setFlag(QGraphicsItem.ItemIsSelectable)

        for name in graph.nodes:
            print(name)

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

        for edge in edges:
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

            self.scene_bottom.addLine(QLineF(x_a*scale, y_a*scale, x_b*scale, y_b*scale))

        for node in nodes:
            x_a,y_a = layout[node + "_left"]
            x_b,y_b = layout[node + "_right"]
            pen = QPen(Qt.green)
            pen.setWidth(3)

            self.scene_bottom.addLine(QLineF(x_a*scale, y_a*scale, x_b*scale, y_b*scale), pen)

    def up(self):
        """ Iterate all selected items in the view, moving them forward. """
        items = self.scene.selectedItems()
        for item in items:
            z = item.zValue()
            item.setZValue(z + 1)

    def down(self):
        """ Iterate all selected items in the view, moving them backward. """
        items = self.scene.selectedItems()
        for item in items:
            z = item.zValue()
            item.setZValue(z - 1)


def main():
    app = QApplication(sys.argv)
    w = Window()
    w.show()
    app.exec()


if __name__ == "__main__":
    main()
