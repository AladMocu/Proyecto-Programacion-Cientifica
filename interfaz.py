from PyQt5 import QtCore, QtGui, QtWidgets
from qt_material import apply_stylesheet
import random
import matplotlib
import numpy as np

matplotlib.use('Qt5Agg')
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from logica import solve


# Canvas de graficas
class MyMplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111)
        self.compute_initial_figure()
        FigureCanvas.__init__(self, self.fig)
        self.setParent(parent)
        FigureCanvas.setSizePolicy(self,
                                   QtWidgets.QSizePolicy.Expanding,
                                   QtWidgets.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
    def compute_initial_figure(self):
        pass


class MyDynamicMplCanvas(MyMplCanvas):
    """A canvas that updates itself every second with a new plot."""
    def __init__(self, *args, **kwargs):
        MyMplCanvas.__init__(self, *args, **kwargs)
        timer = QtCore.QTimer(self)
        timer.timeout.connect(self.update_figure)
        self.axes.set_xlabel('Days')
        self.axes.set_ylabel('Population Radio')
        self.axes.set_title('Method')
        self.cur_sol = solve("odeint/ivp-solve", [0, 0, 0, 0, 0, 0, 0], np.arange(0,150))
        self.cur_meth="odeint/ivp-solve"

    def solve_model(self, method, params, max):
        x_range = np.arange(0, max)
        self.cur_meth = method
        self.cur_sol = solve(method, params, x_range)

    def update_figure(self, variables, max):
        x_range = np.arange(0, max)
        s, e, i, r, p = self.cur_sol
        self.axes.clear()
        if variables[0]:
            self.axes.plot(x_range, s, label="s(t)")
        if variables[1]:
            self.axes.plot(x_range, e, label="e(t)")
        if variables[2]:
            self.axes.plot(x_range, i, label="i(t)")
        if variables[3]:
            self.axes.plot(x_range, r, label="r(t)")
        if variables[4]:
            self.axes.plot(x_range, p, label="p(t)")
        self.axes.set_xlabel('Days')
        self.axes.set_ylabel('Population Radio')
        self.axes.legend()
        self.axes.set_title(self.cur_meth)
        self.draw()


# Dise??o de la interfaz usando qt
# noinspection PyAttributeOutsideInit
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.method = ""
        self.parameteres = [0, 0, 0, 0, 0, 0, 0]

    def setupUi(self, ProyectoFinal):
        ProyectoFinal.setObjectName("ProyectoFinal")
        ProyectoFinal.resize(974, 480)
        self.centralwidget = QtWidgets.QWidget(ProyectoFinal)
        self.centralwidget.setObjectName("centralwidget")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.centralwidget)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.widget_2 = QtWidgets.QWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(3)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.widget_2.sizePolicy().hasHeightForWidth())
        self.widget_2.setSizePolicy(sizePolicy)
        self.widget_2.setObjectName("widget_2")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.widget_2)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setContentsMargins(4, 4, 4, 4)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_2 = QtWidgets.QLabel(self.widget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setObjectName("label_2")
        self.verticalLayout.addWidget(self.label_2)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.formLayout = QtWidgets.QFormLayout()
        self.formLayout.setObjectName("formLayout")
        self.kLabel = QtWidgets.QLabel(self.widget_2)
        self.kLabel.setObjectName("kLabel")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.LabelRole, self.kLabel)
        self.kLineEdit = QtWidgets.QLineEdit(self.widget_2)
        self.kLineEdit.setObjectName("kLineEdit")
        self.formLayout.setWidget(0, QtWidgets.QFormLayout.FieldRole, self.kLineEdit)
        self.aiLabel = QtWidgets.QLabel(self.widget_2)
        self.aiLabel.setObjectName("aiLabel")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.LabelRole, self.aiLabel)
        self.aiLineEdit = QtWidgets.QLineEdit(self.widget_2)
        self.aiLineEdit.setObjectName("aiLineEdit")
        self.formLayout.setWidget(1, QtWidgets.QFormLayout.FieldRole, self.aiLineEdit)
        self.aeLabel = QtWidgets.QLabel(self.widget_2)
        self.aeLabel.setObjectName("aeLabel")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.LabelRole, self.aeLabel)
        self.aeLineEdit = QtWidgets.QLineEdit(self.widget_2)
        self.aeLineEdit.setObjectName("aeLineEdit")
        self.formLayout.setWidget(2, QtWidgets.QFormLayout.FieldRole, self.aeLineEdit)
        self.yLabel = QtWidgets.QLabel(self.widget_2)
        self.yLabel.setObjectName("yLabel")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.LabelRole, self.yLabel)
        self.yLineEdit = QtWidgets.QLineEdit(self.widget_2)
        self.yLineEdit.setObjectName("yLineEdit")
        self.formLayout.setWidget(3, QtWidgets.QFormLayout.FieldRole, self.yLineEdit)
        self.bLabel = QtWidgets.QLabel(self.widget_2)
        self.bLabel.setObjectName("bLabel")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.LabelRole, self.bLabel)
        self.bLineEdit = QtWidgets.QLineEdit(self.widget_2)
        self.bLineEdit.setObjectName("bLineEdit")
        self.formLayout.setWidget(4, QtWidgets.QFormLayout.FieldRole, self.bLineEdit)
        self.pLabel = QtWidgets.QLabel(self.widget_2)
        self.pLabel.setObjectName("pLabel")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.LabelRole, self.pLabel)
        self.pLineEdit = QtWidgets.QLineEdit(self.widget_2)
        self.pLineEdit.setObjectName("pLineEdit")
        self.formLayout.setWidget(5, QtWidgets.QFormLayout.FieldRole, self.pLineEdit)
        self.uLabel = QtWidgets.QLabel(self.widget_2)
        self.uLabel.setObjectName("uLabel")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.LabelRole, self.uLabel)
        self.uLineEdit = QtWidgets.QLineEdit(self.widget_2)
        self.uLineEdit.setObjectName("uLineEdit")
        self.formLayout.setWidget(6, QtWidgets.QFormLayout.FieldRole, self.uLineEdit)
        self.verticalLayout_2.addLayout(self.formLayout)
        self.verticalLayout.addLayout(self.verticalLayout_2)
        self.label = QtWidgets.QLabel(self.widget_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.btn_euler_forward = QtWidgets.QPushButton(self.widget_2)
        self.btn_euler_forward.setObjectName("btn_euler_forward")
        self.verticalLayout.addWidget(self.btn_euler_forward)
        self.btn_euler_backward = QtWidgets.QPushButton(self.widget_2)
        self.btn_euler_backward.setObjectName("btn_euler_backward")
        self.verticalLayout.addWidget(self.btn_euler_backward)
        self.btn_euler_modified = QtWidgets.QPushButton(self.widget_2)
        self.btn_euler_modified.setObjectName("btn_euler_modified")
        self.verticalLayout.addWidget(self.btn_euler_modified)
        self.btn_rk_2 = QtWidgets.QPushButton(self.widget_2)
        self.btn_rk_2.setObjectName("btn_rk_2")
        self.verticalLayout.addWidget(self.btn_rk_2)
        self.btn_rk_4 = QtWidgets.QPushButton(self.widget_2)
        self.btn_rk_4.setObjectName("btn_rk_4")
        self.verticalLayout.addWidget(self.btn_rk_4)
        self.btn_ivp = QtWidgets.QPushButton(self.widget_2)
        self.btn_ivp.setObjectName("btn_ivp")
        self.verticalLayout.addWidget(self.btn_ivp)
        self.gridLayout_2.addLayout(self.verticalLayout, 0, 0, 1, 1)
        self.gridLayout_3.addWidget(self.widget_2, 0, 1, 1, 1)
        self.widget = QtWidgets.QWidget(self.centralwidget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(8)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.widget.sizePolicy().hasHeightForWidth())
        self.widget.setSizePolicy(sizePolicy)
        self.widget.setObjectName("widget")
        self.gridLayout = QtWidgets.QGridLayout(self.widget)
        self.gridLayout.setObjectName("gridLayout")
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setSizeConstraint(QtWidgets.QLayout.SetDefaultConstraint)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.param_s = QtWidgets.QRadioButton(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.param_s.sizePolicy().hasHeightForWidth())
        self.param_s.setSizePolicy(sizePolicy)
        self.param_s.setAutoExclusive(False)
        self.param_s.setObjectName("param_s")
        self.horizontalLayout.addWidget(self.param_s)
        self.param_e = QtWidgets.QRadioButton(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.param_e.sizePolicy().hasHeightForWidth())
        self.param_e.setSizePolicy(sizePolicy)
        self.param_e.setAutoExclusive(False)
        self.param_e.setObjectName("param_e")
        self.horizontalLayout.addWidget(self.param_e)
        self.param_i = QtWidgets.QRadioButton(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.param_i.sizePolicy().hasHeightForWidth())
        self.param_i.setSizePolicy(sizePolicy)
        self.param_i.setAutoExclusive(False)
        self.param_i.setObjectName("param_i")
        self.horizontalLayout.addWidget(self.param_i)
        self.param_r = QtWidgets.QRadioButton(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.param_r.sizePolicy().hasHeightForWidth())
        self.param_r.setSizePolicy(sizePolicy)
        self.param_r.setAutoExclusive(False)
        self.param_r.setObjectName("param_r")
        self.horizontalLayout.addWidget(self.param_r)
        self.param_p = QtWidgets.QRadioButton(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.param_p.sizePolicy().hasHeightForWidth())
        self.param_p.setSizePolicy(sizePolicy)
        self.param_p.setAutoExclusive(False)
        self.param_p.setObjectName("param_p")
        self.horizontalLayout.addWidget(self.param_p)
        self.gridLayout.addLayout(self.horizontalLayout, 1, 0, 1, 1)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.frame = QtWidgets.QFrame(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.frame.sizePolicy().hasHeightForWidth())
        self.frame.setSizePolicy(sizePolicy)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        l = QtWidgets.QVBoxLayout(self.frame)
        self.dc = MyDynamicMplCanvas(self.frame, width=5, height=4, dpi=100)
        l.addWidget(self.dc)
        self.verticalLayout_3.addWidget(self.frame)
        self.gridLayout.addLayout(self.verticalLayout_3, 0, 0, 1, 1)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setContentsMargins(100, -1, 100, -1)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.lineEdit = QtWidgets.QLineEdit(self.widget)

        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit.sizePolicy().hasHeightForWidth())
        self.lineEdit.setSizePolicy(sizePolicy)
        self.lineEdit.setObjectName("lineEdit")
        self.lineEdit.setText("0")
        self.lineEdit.setReadOnly(True)
        self.horizontalLayout_2.addWidget(self.lineEdit)
        self.label_3 = QtWidgets.QLabel(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy)
        self.label_3.setObjectName("label_3")
        self.horizontalLayout_2.addWidget(self.label_3)
        self.lineEdit_2 = QtWidgets.QLineEdit(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_2.sizePolicy().hasHeightForWidth())
        self.lineEdit_2.setSizePolicy(sizePolicy)
        self.lineEdit_2.setObjectName("lineEdit_2")
        self.horizontalLayout_2.addWidget(self.lineEdit_2)
        self.label_4 = QtWidgets.QLabel(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_4.sizePolicy().hasHeightForWidth())
        self.label_4.setSizePolicy(sizePolicy)
        self.label_4.setObjectName("label_4")
        self.horizontalLayout_2.addWidget(self.label_4)
        self.gridLayout.addLayout(self.horizontalLayout_2, 4, 0, 1, 1)
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_3.setSizeConstraint(QtWidgets.QLayout.SetMinimumSize)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.label_5 = QtWidgets.QLabel(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_5.sizePolicy().hasHeightForWidth())
        self.label_5.setSizePolicy(sizePolicy)
        self.label_5.setObjectName("label_5")
        self.horizontalLayout_3.addWidget(self.label_5)
        self.sim_time = QtWidgets.QLCDNumber(self.widget)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.sim_time.sizePolicy().hasHeightForWidth())
        self.sim_time.setSizePolicy(sizePolicy)
        self.sim_time.setSmallDecimalPoint(False)
        self.sim_time.setObjectName("sim_time")
        self.horizontalLayout_3.addWidget(self.sim_time)
        self.gridLayout.addLayout(self.horizontalLayout_3, 3, 0, 1, 1)
        self.gridLayout_3.addWidget(self.widget, 0, 0, 1, 1)
        ProyectoFinal.setCentralWidget(self.centralwidget)
        self.toolbar = ProyectoFinal.addToolBar('Tools')
        self.exitAction = QtWidgets.QAction(QtGui.QIcon('./assets/exit.png'), 'Exit', self)
        self.exitAction.triggered.connect(ProyectoFinal.close)
        self.toolbar.addAction(self.exitAction)

        self.importAction = QtWidgets.QAction(QtGui.QIcon('./assets/import.png'), 'Import', self)
        self.importAction.triggered.connect(lambda: self.import_data())
        self.toolbar.addAction(self.importAction)

        self.exportAction = QtWidgets.QAction(QtGui.QIcon('./assets/export.png'), 'Export', self)
        self.exportAction.triggered.connect(lambda: self.export_data())
        self.toolbar.addAction(self.exportAction)

        self.statusbar = QtWidgets.QStatusBar(ProyectoFinal)
        self.statusbar.setObjectName("statusbar")
        ProyectoFinal.setStatusBar(self.statusbar)
        self.actionSalir = QtWidgets.QAction(ProyectoFinal)
        self.actionSalir.setObjectName("actionSalir")
        self.actionImportar = QtWidgets.QAction(ProyectoFinal)
        self.actionImportar.setObjectName("actionImportar")
        self.actionExportar = QtWidgets.QAction(ProyectoFinal)
        self.actionExportar.setObjectName("actionExportar")
        # self.menubar.addAction(self.menuStyles.menuAction())
        # self.menubar.addAction(self.menuImportar.menuAction())
        # self.menubar.addAction(self.menuExportar.menuAction())

        self.retranslateUi(ProyectoFinal)
        QtCore.QMetaObject.connectSlotsByName(ProyectoFinal)

    def retranslateUi(self, ProyectoFinal):
        _translate = QtCore.QCoreApplication.translate
        ProyectoFinal.setWindowTitle(_translate("ProyectoFinal", "ProyectoFinal"))
        self.label_2.setText(_translate("ProyectoFinal", "Parametros"))
        self.kLabel.setText(_translate("ProyectoFinal", "??"))
        self.aiLabel.setText(_translate("ProyectoFinal", "??????"))
        self.aeLabel.setText(_translate("ProyectoFinal", "??????"))
        self.yLabel.setText(_translate("ProyectoFinal", "??"))
        self.bLabel.setText(_translate("ProyectoFinal", "??"))
        self.pLabel.setText(_translate("ProyectoFinal", "??"))
        self.uLabel.setText(_translate("ProyectoFinal", "??"))
        self.label.setText(_translate("ProyectoFinal", "Metodos de solucion"))
        self.btn_euler_forward.setText(_translate("ProyectoFinal", "Euler Adelante"))
        self.btn_euler_backward.setText(_translate("ProyectoFinal", "Euler Atras"))
        self.btn_euler_modified.setText(_translate("ProyectoFinal", "Euler Modificado"))
        self.btn_rk_2.setText(_translate("ProyectoFinal", "Runge-Kutta 2"))
        self.btn_rk_4.setText(_translate("ProyectoFinal", "Runge-Kutta 4"))
        self.btn_ivp.setText(_translate("ProyectoFinal", "odeint/ivp-solve"))
        self.param_s.setText(_translate("ProyectoFinal", "s(t)"))
        self.param_e.setText(_translate("ProyectoFinal", "e(t)"))
        self.param_i.setText(_translate("ProyectoFinal", "i(t)"))
        self.param_r.setText(_translate("ProyectoFinal", "r(t)"))
        self.param_p.setText(_translate("ProyectoFinal", "p(t)"))
        self.label_3.setText(_translate("ProyectoFinal", "-"))
        self.label_4.setText(_translate("ProyectoFinal", "Dias"))
        self.label_5.setText(_translate("ProyectoFinal", "Duracion Simulacion:"))
        # self.menuStyles.setTitle(_translate("ProyectoFinal", "Salir"))
        # self.menuImportar.setTitle(_translate("ProyectoFinal", "Importar"))
        # self.menuExportar.setTitle(_translate("ProyectoFinal", "Exportar"))
        self.actionSalir.setText(_translate("ProyectoFinal", "Salir"))
        self.actionImportar.setText(_translate("ProyectoFinal", "Importar"))
        self.actionImportar.setToolTip(_translate("ProyectoFinal", "Importar"))
        self.actionExportar.setText(_translate("ProyectoFinal", "Exportar"))
        self.actionExportar.setToolTip(_translate("ProyectoFinal", "Exportar"))
        self.kLineEdit.setPlaceholderText("0.050")
        self.aiLineEdit.setPlaceholderText("0.005")
        self.aeLineEdit.setPlaceholderText("0.650")
        self.yLineEdit.setPlaceholderText("0.000")
        self.bLineEdit.setPlaceholderText("0.100")
        self.pLineEdit.setPlaceholderText("0.080")
        self.uLineEdit.setPlaceholderText("0.020")
        self.lineEdit_2.setPlaceholderText("150")
        self.add_actions()

    def add_actions(self):
        # acciones
        # self.menuStyles.clicked.connect(lambda: self.salir())
        # self.menuExportar.clicked.connect(lambda: self.export_data())
        # self.menuExportar.clicked.connect(lambda: self.export_data())
        # variables
        self.param_s.toggled.connect(lambda: self.update_plot(False))
        self.param_e.toggled.connect(lambda: self.update_plot(False))
        self.param_i.toggled.connect(lambda: self.update_plot(False))
        self.param_r.toggled.connect(lambda: self.update_plot(False))
        self.param_p.toggled.connect(lambda: self.update_plot(False))
        # metodos
        self.btn_euler_forward.clicked.connect(lambda: self.change_method("Euler Forward"))
        self.btn_euler_backward.clicked.connect(lambda: self.change_method("Euler Backward"))
        self.btn_euler_modified.clicked.connect(lambda: self.change_method("Euler Modified"))
        self.btn_rk_2.clicked.connect(lambda: self.change_method("Runge-Kutta 2"))
        self.btn_rk_4.clicked.connect(lambda: self.change_method("Runge-Kutta 4"))
        self.btn_ivp.clicked.connect(lambda: self.change_method("odeint/ivp-solve"))
        # parametros
        self.kLineEdit.setValidator(QtGui.QDoubleValidator(-1000.0, 1000.0, 2))
        self.kLineEdit.setValidator(QtGui.QDoubleValidator(-1000.0, 1000.0, 2))
        self.aiLineEdit.setValidator(QtGui.QDoubleValidator(-1000.0, 1000.0, 2))
        self.aeLineEdit.setValidator(QtGui.QDoubleValidator(-1000.0, 1000.0, 2))
        self.yLineEdit.setValidator(QtGui.QDoubleValidator(-1000.0, 1000.0, 2))
        self.bLineEdit.setValidator(QtGui.QDoubleValidator(-1000.0, 1000.0, 2))
        self.pLineEdit.setValidator(QtGui.QDoubleValidator(-1000.0, 1000.0, 2))
        self.uLineEdit.setValidator(QtGui.QDoubleValidator(-1000.0, 1000.0, 2))


    def change_method(self, newMethod):
        self.method = newMethod
        self.update_plot(True)

    def update_plot(self, recalculate):
        t_vars = [self.param_s.isChecked(), self.param_e.isChecked(), self.param_i.isChecked(),
                  self.param_r.isChecked(), self.param_p.isChecked()]
        self.parameteres = [
            float(self.kLineEdit.text() or "0.050"),
            float(self.aiLineEdit.text() or "0.005"),
            float(self.aeLineEdit.text() or "0.650"),
            float(self.yLineEdit.text() or "0.000"),
            float(self.bLineEdit.text() or "0.100"),
            float(self.pLineEdit.text() or "0.080"),
            float(self.uLineEdit.text() or "0.020")]

        print("Method: ", self.method)
        print("Parameters: ", self.parameteres)
        print("Variables: ", t_vars)
        if recalculate:
            self.dc.solve_model(self.method,self.parameteres,float(self.lineEdit_2.text() or "150"))
        self.dc.update_figure(t_vars, float(self.lineEdit_2.text() or "150"))

    def import_data(self):
        path = QtWidgets.QFileDialog.getOpenFileName(self, 'Abrir un archivo', '', 'Poputation Simulation files (*.sd)')
        if path != ('', ''):
            print("File path : " + path[0])
            self.parameteres = np.fromfile(path[0], dtype=np.float32)
            self.kLineEdit.setText(str(self.parameteres[0]))
            self.aiLineEdit.setText(str(self.parameteres[1]))
            self.aeLineEdit.setText(str(self.parameteres[2]))
            self.yLineEdit.setText(str(self.parameteres[3]))
            self.bLineEdit.setText(str(self.parameteres[4]))
            self.pLineEdit.setText(str(self.parameteres[5]))
            self.uLineEdit.setText(str(self.parameteres[6]))

    def export_data(self):
        self.parameteres = [
            float(self.kLineEdit.text() or "0"),
            float(self.aiLineEdit.text() or "0"),
            float(self.aeLineEdit.text() or "0"),
            float(self.yLineEdit.text() or "0"),
            float(self.bLineEdit.text() or "0"),
            float(self.pLineEdit.text() or "0"),
            float(self.uLineEdit.text() or "0")]

        name = QtWidgets.QFileDialog.getSaveFileName(self, 'Exportar datos', '', 'Poputation Simulation files (*.sd)')
        if name != ('', ''):
            print(name[0])
            data = np.array(self.parameteres, dtype=np.float32)
            data.tofile(name[0])

    def update_duracion(self):
        self.sim_time.value = 32


def draw():
    import sys
    app = QtWidgets.QApplication(sys.argv)
    ProyectoFinal = QtWidgets.QMainWindow()
    apply_stylesheet(app, theme='dark_red.xml')
    ui = MainWindow()
    ui.setupUi(ProyectoFinal)
    ProyectoFinal.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    draw()
