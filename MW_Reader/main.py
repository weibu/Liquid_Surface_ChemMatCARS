import sys

# import PyQt4 QtCore and QtGui modules
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from mainwindow import MainWindow

if __name__ == '__main__':

    # create application
    app = QApplication(sys.argv)
    app.setApplicationName('MW-XReader')
    

    # create widget
    w = MainWindow()
    w.setWindowTitle('MW-XReader')
    w.setWindowIcon(QIcon('logo.png'))
    w.show()

    # connection
    QObject.connect(app, SIGNAL('lastWindowClosed()'), app, SLOT('quit()'))

    # execute application
    sys.exit(app.exec_())
