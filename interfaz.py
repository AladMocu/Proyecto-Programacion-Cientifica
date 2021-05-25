from PySide6.QtWidgets import QApplication, QMainWindow
from PySide6.QtUiTools import QUiLoader
from PySide6.QtGui import QFontDatabase
from qt_material import QtStyleTools


########################################################################
class RuntimeStylesheets(QMainWindow, QtStyleTools):
    # ----------------------------------------------------------------------
    def __init__(self):
        """"""
        super().__init__()
        self.main = QUiLoader().load('window.ui', self)
        self.apply_stylesheet(self.main, theme='dark_teal.xml')

def draw():
    app = QApplication()
    frame = RuntimeStylesheets()
    frame.main.show()
    app.exec_()



