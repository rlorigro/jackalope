import os.path

import sys

from PyQt5.QtCore import QUrl
from PyQt5.QtWebEngineWidgets import QWebEngineView
from PyQt5.QtWidgets import (
    QApplication,
    QGraphicsScene,
    QHBoxLayout,
    QWidget,
)

#
# Adapted from https://www.pythonguis.com/tutorials/pyqt-qgraphics-vector-graphics/
#

sys.argv.append("--disable-web-security")

class Window(QWidget):
    def __init__(self):
        super().__init__()

        # Defining a scene rect of 400x200, with its origin at 0,0.
        # If we don't set this on creation, we can set it later with .setSceneRect
        self.project_dir = os.path.dirname(os.path.dirname(__file__))
        self.igv_html = os.path.join(self.project_dir, "html/igv.html")

        print(self.igv_html)

        self.view_top = QWebEngineView()
        self.view_top.resize(800,400)
        # self.view_top.load(QUrl("http://www.google.com/"))
        self.view_top.load(QUrl.fromLocalFile(self.igv_html))
        self.view_top.show()
        self.view_top.loadFinished.connect(self.load_finished)

        hbox = QHBoxLayout(self)
        hbox.addWidget(self.view_top)

        self.setLayout(hbox)

    def load_finished(self):
        print("done")


app = QApplication(sys.argv)

w = Window()
w.resize(800,400)
w.show()

app.exec()