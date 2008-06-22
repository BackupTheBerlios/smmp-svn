#!/usr/bin/env python
# gui_example.py
#
# Copyright 2007       Frank Eisenmenger, U.H.E. Hansmann,
#                      Jan H. Meinke, Sandipan Mohanty
#
# Display the Calpha trace and energy of a Monte Carlo simulation of a protein using Qt4.
import sys
sys.path.append("../../")

import universe
import protein
import algorithms
from PyQt4 import  QtCore, QtGui

class CAlphaTraceWidget(QtGui.QWidget):
    """Displays a 2D projection of the Calpha trace of the protein p. The 
    projection simply drops the z coordinate.
    """
    def __init__(self, protein, parent = 0):
        """Initialize the widget.
        
        @param p the protein.Protein to be displayed
        @param parent parent widget if applicable.
        """
        QtGui.QWidget.__init__(self)
        self.protein = protein
        self.calphaList = []
        for a in self.protein.atoms():
            if a.name().lower().strip() == 'ca':
                self.calphaList.append(a)
        self.centralAtom = self.calphaList[len(self.calphaList) / 2]
         
    def paintEvent(self, event):
        """Draw little circles in the position of the Calpha atoms of p."""
        diameter = 2
        avgCalphaDistance = 4
        zmin = -avgCalphaDistance * len(self.protein)
        zmax = avgCalphaDistance * len(self.protein)
        
        self.calphaList.sort(self.atomZSort)

        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.Antialiasing, 1)
        painter.translate(self.width() / 2, self.height() / 2)
        xscale = 1.5 * self.width() / (avgCalphaDistance * len(self.protein))
        yscale = 1.5 * self.height() / (avgCalphaDistance * len(self.protein))
        scale = min(xscale, yscale)
        painter.scale(scale, scale) 
        myGradient = QtGui.QRadialGradient(QtCore.QPointF(0, 0), 72)
        myGradient.setColorAt(0, QtCore.Qt.white)
        myGradient.setColorAt(1, QtCore.Qt.black)
        painter.setBrush(QtGui.QBrush(myGradient))

        painter.translate(-self.centralAtom.position()[0], -self.centralAtom.position()[1])

        
        for a in self.calphaList:
            pos = a.position()
            if pos[2] < zmin:
                zscale = 0.1
            elif pos[2] > zmax:
                zscale = 5.0
            else:
                zscale = (pos[2] - zmin) / (zmax - zmin) * 4.9 + 0.1

            painter.save()
            painter.translate(pos[0], pos[1])
            painter.scale(zscale * diameter / 100, zscale * diameter / 100)
            painter.drawEllipse(-50, -50, 100, 100)
            painter.restore()
            
        painter.setWorldMatrixEnabled(0)
        painter.drawText(5,self.height() -5 , QtCore.QString("E=%s" % self.protein.energy()))
        painter.setWorldMatrixEnabled(1)
        
    def atomZSort(self, a, b):
        za = a.position()[2]
        zb = b.position()[2]
        if za < zb:
            return -1
        elif za > zb:
            return 1
        else:
            return 0


class WorkThread(QtCore.QThread):
    def __init__(self, f):
        QtCore.QThread.__init__(self)
        self.function= f
        
    def run(self):
        import os.path
        from os import system
        for i in range(0, 1000):
            print "Performing sweep %s." % i
            self.function()
            try:
                print os.path.getmtime('best.png'), os.path.getmtime('best.pdb')
                if os.path.getmtime('best.png') < os.path.getmtime('best.pdb'):
                    os.system('pymol -c best.pml')
            except os.error,e:
                os.system('pymol -c best.pml')

            self.emit(QtCore.SIGNAL('finishedSweep(int)'), i)
    
class TemperatureDial(QtGui.QDial):

    def __init__(self, Tmin, Tmax):
        QtGui.QDial.__init__(self)
        self.setMinimum(Tmin)
        self.setMaximum(Tmax)
        self.setNotchTarget(60)
        self.setNotchesVisible(1)
        
class MainWidget(QtGui.QWidget):

        def __init__(self):
            QtGui.QWidget.__init__(self)
            self.resize(640,480)
            self.setWindowTitle("MC Simulation of Protein A (1BDD)")
            topLevelLayout = QtGui.QHBoxLayout()
            sideBarLayout = QtGui.QVBoxLayout()

            self.maxSweeps = 1000
            self.myU = universe.Universe()
            self.p = protein.Protein('../1bdd.seq')
            self.myMC = algorithms.CanonicalMonteCarlo(self.myU, 1, self.maxSweeps)
            self.myWorkThread = WorkThread(self.myMC.sweep)

            self.labelList = [['E', self.p.energy], ['R<sub>gyr</sub>', self.p.rgyr, 0], ['hb' ,self.p.hbond]]

            myCAlphaTraceWidget = CAlphaTraceWidget(self.p)
            topLevelLayout.addWidget(myCAlphaTraceWidget, 1)

            self.createSideBar(sideBarLayout)
            topLevelLayout.addLayout(sideBarLayout)

            self.setLayout(topLevelLayout)
            
            self.myWorkThread.start()

            self.connect(self.myWorkThread, QtCore.SIGNAL("finishedSweep(int)"), self.finishedSweep)
            self.connect(self.myWorkThread, QtCore.SIGNAL("finishedSweep(int)"), self.myProgressBar, QtCore.SLOT("setValue(int)"))

        def createSideBar(self, sideBarLayout):
        
            myTDial = TemperatureDial(1, 1200)
            myTDial.setValue(self.myU.temperature())
            sideBarLayout.addWidget(myTDial)
            temperatureLabel = QtGui.QLabel("%s" % myTDial.value())
            self.connect(myTDial, QtCore.SIGNAL("valueChanged(int)"), temperatureLabel, QtCore.SLOT("setNum(int)"))
            self.connect(myTDial, QtCore.SIGNAL("valueChanged(int)"), self.myU.setTemperature)
            sideBarLayout.addWidget(temperatureLabel)
            sideBarLayout.addStretch(1)
            self.labels = {}
            for name in self.labelList:
                if len(name) > 2:
                    labelText = "%s: %8.3f" % (name[0], name[1]()[name[2]])
                else:
                    labelText = "%s: %8.3f" % (name[0], name[1]())
                self.labels[name[0]] = QtGui.QLabel(labelText)
                sideBarLayout.addWidget(self.labels[name[0]])
            sideBarLayout.addStretch(2)
            self.bestImageLabel = QtGui.QLabel()
            self.bestImageLabel.setText('Waiting...')
            self.bestImageLabel.setToolTip('Waiting ...')
            sideBarLayout.addWidget(self.bestImageLabel)

            self.myProgressBar = QtGui.QProgressBar()
            self.myProgressBar.setRange(0, self.maxSweeps - 1)
            self.myProgressBar.setValue(0)
            
            sideBarLayout.addWidget(self.myProgressBar)
            
        def finishedSweep(self, sweep):
            for name in self.labelList:
                if len(name) > 2:
                    labelText = "%s: %8.3f" % (name[0], name[1]()[name[2]])
                else:
                    labelText = "%s: %8.3f" % (name[0], name[1]())
                self.labels[name[0]].setText(labelText)
            bestPix = QtGui.QPixmap("best.png")
            scaledPix = bestPix.scaledToWidth(100)
            self.bestImageLabel.setPixmap(scaledPix)
            self.bestImageLabel.setToolTip('<img src="best.png"> E=%s at t=%s' % self.myMC.minimumEnergy())
            self.setWindowIcon(QtGui.QIcon(scaledPix))
            self.update()
                
if __name__ == "__main__":
    myApp = QtGui.QApplication(sys.argv)
    myMain = MainWidget()
    myMain.show()
    sys.exit(myApp.exec_())
