import sys
from PySide6 import QtWidgets
from qt_material import apply_stylesheet

#primer archivo
def draw():
    app = QtWidgets.QApplication(sys.argv)
    window = QtWidgets.QMainWindow()
    window.resize(700, 700)
    window.show()
    app.exec_()


if __name__ == "__main__":
    print("por favor ejecuta main.py")
