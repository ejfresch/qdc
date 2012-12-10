#!/usr/bin/python

# Copyright 2009-2012 - Luca Freschi <l.freschi@gmail.com>                                                                       
# This file is part of QDC.                                                                                                  
# QDC is free software: you can redistribute it and/or modify                                                               
# it under the terms of the GNU General Public License as published by                                                      
# the Free Software Foundation, either version 3 of the License, or                                                         
# (at your option) any later version.  

# This program is distributed in the hope that it will be useful,                                                           
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                                            
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                             
# GNU General Public License for more details.                                                                              

# You should have received a copy of the GNU General Public License                                                         
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 

import sys
import commands
import os
from os.path import isfile
import csv
import re

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib import rc


from PyQt4 import QtGui
from PyQt4 import QtCore

import import_sbml as i_sbml
import export_sbml as e_sbml


class Trd(QtCore.QThread):
	def __init__(self, filename, num_of_simulations, parent=None):
		QtCore.QThread.__init__(self, parent)
		
		self.f=filename
		self.n=num_of_simulations
		
		
		
	def run(self):
		
		cmd='python taxi_driver.py '+str(self.f)+' '+str(self.n)+' > qdc.out'
		status, out = commands.getstatusoutput(cmd)
	

		
		

		
	
class Plot(QtCore.QThread):
	
	def __init__(self, csv_variable, dic, graphs, canvas, axes, labels, parent=None):

	# grey scale 
	#def __init__(self, csv_variable, dic, graphs, canvas, axes, labels, count_dash, dashes, parent=None):
		QtCore.QThread.__init__(self, parent)
		
		self.csv_variable=csv_variable
		self.dic=dic
		self.graphs=graphs
		self.canvas=canvas
		self.axes=axes
		self.labels=labels
		#grey scalep		
		#self.count_dash=count_dash
		#self.dashes=dashes
		
	def run(self):

		self.query=str(self.csv_variable)	
		if self.query != 'Choose variable...':
			self.array=[]
			
			self.labels.append(self.query)
			
			for row in self.dic:
				self.array.append(float(row[self.query]))
			
			length_array=len(self.array)
			self.x=range(length_array)
				

			#grey scale
			#ind=(self.count_dash % len(self.dashes))
			#l1=self.axes.plot(self.x, self.array, linestyle=self.dashes[ind], color='black')
			
			l1=self.axes.plot(self.x, self.array)
			self.graphs.append(l1)
			self.axes.legend(self.graphs,self.labels)
			
			
			
			
		
			
	





class SimGui(QtGui.QWidget):
	def __init__(self, parent=None):
		QtGui.QWidget.__init__(self, parent)
		
		global gt1
		t1=os.getpid()
		#variables		
		self.filename=False		
		self.csv_filename=False
		self.csv_variable=False
		#self.open_flag=1
		self.dic=[]
		self.query=False
		self.file_csv=False
		self.names=[]
		self.graphs=[]
		self.labels=[]
		rc("font", size=10)
		
		
		#GUI
		self.setWindowTitle('QDC')
		self.setWindowIcon(QtGui.QIcon('icons/lab.png'))
		self.file_inspector=QtGui.QLabel('No file loaded', self)
		self.load_button=QtGui.QPushButton("load", self)
		self.load_button.setIcon(QtGui.QIcon('icons/open.png'))		
		self.text_box=QtGui.QTextEdit()
		self.save_button=QtGui.QPushButton("save", self)
		self.save_button.setIcon(QtGui.QIcon('icons/save.png'))
		self.save_as_button=QtGui.QPushButton("save as", self)
		self.save_as_button.setIcon(QtGui.QIcon('icons/save_as.png'))		
		self.import_sbml_button=QtGui.QPushButton("import SBML", self)
		self.import_sbml_button.setIcon(QtGui.QIcon('icons/import.png'))
		self.export_sbml_button=QtGui.QPushButton("export SBML", self)	
		self.export_sbml_button.setIcon(QtGui.QIcon('icons/export.png'))	
		self.quit_button=QtGui.QPushButton("quit", self)
		self.quit_button.setIcon(QtGui.QIcon('icons/close.png'))	
		self.log_inspector=QtGui.QTextEdit()	
		self.start_button=QtGui.QPushButton("Start simulations", self)
		self.start_button.setIcon(QtGui.QIcon('icons/run.png'))	
		self.stop_button=QtGui.QPushButton("Stop", self)
                self.stop_button.setIcon(QtGui.QIcon('icons/del.png'))
		self.status_inspector=QtGui.QLabel('Ready', self)
		self.loadcsv_button=QtGui.QPushButton("load csv", self)
		self.loadcsv_button.setIcon(QtGui.QIcon('icons/open.png'))		
		self.csv_file=QtGui.QLabel('No file loaded',self)
		self.csv_variables=QtGui.QComboBox()
		self.addplot_button=QtGui.QPushButton("add plot", self)
		self.addplot_button.setIcon(QtGui.QIcon('icons/add.png'))
		self.clearplot_button=QtGui.QPushButton("clear plot", self)
		self.clearplot_button.setIcon(QtGui.QIcon('icons/clear.png'))
		self.saveplot_button=QtGui.QPushButton("save plot", self)
		self.saveplot_button.setIcon(QtGui.QIcon('icons/save_as.png'))
		self.sim_label_1=QtGui.QLabel("simulations", self)
		self.n_of_simulations=QtGui.QSpinBox(self)
		self.n_of_simulations.setValue(1)
		self.n_of_simulations.setMinimum(1)
		self.n_of_simulations.setRange(1, 100000)
		

		self.qwidget_figure = QtGui.QWidget()
        
        	
        	self.dpi = 100
        	self.fig = Figure((4.5, 2.5), dpi=self.dpi)
		
        	self.canvas = FigureCanvas(self.fig)
        	self.canvas.setParent(self.qwidget_figure)
		self.axes = self.fig.add_subplot(111)



		self.tabs=QtGui.QTabWidget()
		
		self.qwidget_edit=QtGui.QWidget()		
		tab_1=QtGui.QGridLayout()
		tab_1.setSpacing(5)
		tab_1.addWidget(self.file_inspector, 1, 0, 1, 4)
		tab_1.addWidget(self.load_button, 1, 4, 1, 1)
		tab_1.addWidget(self.text_box, 2, 0, 5, 5)
		tab_1.addWidget(self.save_button, 7, 0, 1, 1)
		tab_1.addWidget(self.save_as_button, 7, 1, 1, 1)
		tab_1.addWidget(self.import_sbml_button, 7, 2, 1, 1)
		tab_1.addWidget(self.export_sbml_button, 7, 3, 1, 1)
		tab_1.addWidget(self.quit_button, 7, 4, 1, 1)
		self.qwidget_edit.setLayout(tab_1)
		self.tabs.addTab(self.qwidget_edit, "edit")
		
		
		self.qwidget_simulation=QtGui.QWidget()		
		tab_2=QtGui.QGridLayout()
		tab_2.setSpacing(5)
		tab_2.addWidget(self.start_button, 2, 2, 1, 3)
		tab_2.addWidget(self.stop_button, 2, 0, 1, 2)
		tab_2.addWidget(self.status_inspector, 1, 4, 1 ,1)
		tab_2.addWidget(self.sim_label_1, 1, 0, 1, 1)
		tab_2.addWidget(self.n_of_simulations, 1, 1, 1, 1)


		tab_2.addWidget(self.log_inspector, 3, 0, 3, 5)
		self.qwidget_simulation.setLayout(tab_2)
		self.tabs.addTab(self.qwidget_simulation, "simulation")

		self.qwidget_plot=QtGui.QWidget()		
		tab_3=QtGui.QGridLayout()
		tab_3.setSpacing(5)
		tab_3.addWidget(self.csv_file, 1, 0, 1, 3)
		tab_3.addWidget(self.loadcsv_button, 1, 3, 1, 1)
		tab_3.addWidget(self.csv_variables, 2, 0, 1, 3)
		tab_3.addWidget(self.addplot_button, 2, 3, 1, 1)
		tab_3.addWidget(self.clearplot_button, 7, 0, 1, 1)
		tab_3.addWidget(self.saveplot_button, 7, 1, 1, 1)
		tab_3.addWidget(self.qwidget_figure, 3, 0, 3, 4)
		self.qwidget_plot.setLayout(tab_3)
		self.tabs.addTab(self.qwidget_plot, "plot")
			
		
		self.layout =QtGui.QVBoxLayout()
		self.layout.addWidget(self.tabs)
		self.setLayout(self.layout)

		self.resize(650, 450)
		self.center()		
	
		#connections
		self.connect(self.load_button, QtCore.SIGNAL('clicked()'), self.open_file)
		self.connect(self.save_button, QtCore.SIGNAL('clicked()'), self.save_file)
		self.connect(self.save_as_button, QtCore.SIGNAL('clicked()'), self.save_file_as)
		self.connect(self.quit_button, QtCore.SIGNAL('clicked()'), self.close_app)
		self.connect(self.start_button, QtCore.SIGNAL('clicked()'), self.simulate)
		self.connect(self.csv_file, QtCore.SIGNAL("currentIndexChanged(int)"),self.onIndexChanged)		
		self.connect(self.loadcsv_button,QtCore.SIGNAL('clicked()'), self.load_csv)		
		self.connect(self.addplot_button,QtCore.SIGNAL('clicked()'), self.addplot)
		self.connect(self.csv_variables, QtCore.SIGNAL("currentIndexChanged(int)"),self.onVariableChanged)
		self.connect(self, QtCore.SIGNAL('closeEmitApp()'), QtCore.SLOT('close()') )
		self.connect(self.clearplot_button, QtCore.SIGNAL('clicked()'), self.clear_plot)
		self.connect(self.saveplot_button, QtCore.SIGNAL('clicked()'), self.save_plot)
		self.connect(self.import_sbml_button, QtCore.SIGNAL('clicked()'), self.import_sbml)
		self.connect(self.export_sbml_button, QtCore.SIGNAL('clicked()'), self.export_sbml)
		self.connect(self.stop_button, QtCore.SIGNAL('clicked()'), self.stop_sim)

		
		

	#functions

	
		


	def center(self):
		screen=QtGui.QDesktopWidget().screenGeometry()
		size=self.geometry()
		self.move((screen.width()-size.width())/2,(screen.height()-size.height())/2)


	def open_file(self):
		c_file=QtGui.QFileDialog.getOpenFileName(self, 'Open file', '.')
		if c_file and isfile(c_file):
			inp=open(c_file)
			contents=inp.read()
			inp.close()
			self.text_box.setText(contents)
			self.file_inspector.setText(c_file)
			self.filename=c_file
			

	def save_file(self):
		if str(self.filename) == 'False':
			n_file=QtGui.QFileDialog(self)
			d_file=n_file.getSaveFileName()
			if d_file:
				n_out=open(d_file,'w')
				my_text=str(self.text_box.toPlainText())						
				n_out.write(my_text)
				n_out.close()
			self.file_inspector.setText(d_file)
			self.filename=d_file

		elif self.filename and isfile(self.filename):
			out=open(self.filename,'w')
			my_text=str(self.text_box.toPlainText())		
										
			out.write(my_text)
			out.close()
		

	def save_file_as(self):
		n_file=QtGui.QFileDialog(self)
		d_file=n_file.getSaveFileName()
		if d_file:
			n_out=open(d_file,'w')
			n_out.write(self.text_box.toPlainText())
			n_out.close()


	def simulate(self):
		if str(self.filename)=='False':
			self.log_inspector.setText(':: No file loaded!')
			return
		self.start_button.setEnabled(False)
		self.status_inspector.setText('Working...')
		self.log_inspector.setText('')
		QtGui.QApplication.processEvents()
		num_of_simulations=self.n_of_simulations.value()	

		
		self.thread=Trd(self.filename, self.n_of_simulations.value())
		self.thread.start()
		
		self.connect(self.thread, QtCore.SIGNAL("finished()"), self.updateUi)
		
		
		return
	


	def stop_sim(self):
		os.popen("killall "+"./parser > /dev/null 2>&1")
		os.popen("killall "+"./engine > /dev/null 2>&1")
		self.log_inspector.setText('')
		self.status_inspector.setText('Ready')
                self.start_button.setEnabled(True)




	def updateUi(self):
		self.log_inspector.setText('')
		ctrl_file="qdc.out"
		if ctrl_file and isfile(ctrl_file):
			f=open('qdc.out', 'r')
			for line in f:

				self.log_inspector.append(line.strip())
			f.close()
		self.status_inspector.setText('Ready')	
		self.start_button.setEnabled(True)
	

	def addplot(self):
		self.thread2=Plot(self.csv_variable, self.dic, self.graphs, self.canvas, self.axes, self.labels)
		# grey scale
		#self.thread2=Plot(self.csv_variable, self.dic, self.graphs, self.canvas, self.axes, self.labels, self.count_dash,self.dashes)
		self.thread2.start()
		self.connect(self.thread2, QtCore.SIGNAL("finished()"), self.updateCanvas)
		#grey scale
		#self.count_dash=self.count_dash+1		
		self.clearplot_button.setEnabled(False)
		self.addplot_button.setEnabled(False)
		self.loadcsv_button.setEnabled(False)
		self.saveplot_button.setEnabled(False)
		

	def updateCanvas(self):
		self.canvas.draw()
		self.clearplot_button.setEnabled(True)
		self.addplot_button.setEnabled(True)
		self.loadcsv_button.setEnabled(True)
		self.saveplot_button.setEnabled(True)
	
	

	def onIndexChanged(self, index):
     		self.csv_filename=self.csv_file.currentText()




	def load_csv(self):
		self.csv_filename=QtGui.QFileDialog.getOpenFileName(self, 'Open file', '.')
				
		if self.csv_filename and isfile(self.csv_filename):
			self.csv_variables.clear()
			self.file_csv=open(self.csv_filename, 'r')
			self.names=csv.reader(self.file_csv).next()
			self.csv_variables.addItem('Choose variable...')
			for elem in self.names:
				self.csv_variables.addItem(elem)
			self.dic=csv.DictReader(self.file_csv, fieldnames=self.names)
			self.csv_file.setText(str(self.csv_filename))
			#grey scale
			#self.dashes=['-',':','--','-.']
			#self.count_dash=0;
			
			return self.file_csv, self.dic, self.names
	
	def onVariableChanged(self, index):
     		self.csv_variable=self.csv_variables.currentText()
		self.file_csv.close()
		if self.csv_filename and isfile(self.csv_filename):
			self.file_csv=open(self.csv_filename, 'r')
			self.names=csv.reader(self.file_csv).next()
			self.dic=csv.DictReader(self.file_csv, fieldnames=self.names)
			return self.dic

	
		

	def close_app(self):
		self.csv_file.close()

		self.close()



	def clear_plot(self):
		self.axes.clear()
		self.labels=[]
		self.graphs=[]
		self.canvas.draw()
		self.csv_variables.setCurrentIndex(0)
		QtGui.QApplication.processEvents()



	def save_plot(self):
		#file_choices = "PNG (*.png)|*.png"
		file_choices = "*.png"
		
		img = unicode(QtGui.QFileDialog.getSaveFileName(self, 'Save file', '.', file_choices))
		if img:
			self.fig.savefig(img, dpi=self.dpi)



	def import_sbml(self):
		sbml_i_file=QtGui.QFileDialog.getOpenFileName(self, 'Open file', '.')
		if sbml_i_file and isfile(sbml_i_file):
			contents=i_sbml.import_sbml(str(sbml_i_file))
			self.text_box.setText(contents)
			self.filename=False
			self.file_inspector.setText("No file loaded")		
			

	def export_sbml(self):
				
		sbml_e_file_choose=QtGui.QFileDialog(self)
		sbml_e_file=sbml_e_file_choose.getSaveFileName()
		to_parse=self.text_box.toPlainText()
		
		my_to_parse=str(to_parse)
		
		pat=re.compile('\n')
		list_to_parse=re.split(pat, my_to_parse)
		
		my_out_sbml=e_sbml.export_sbml(list_to_parse)
		
		if sbml_e_file:
			f_out=open(sbml_e_file,'w')
			f_out.write(str(my_out_sbml))
			f_out.close()


		
		
def main(self):	
	application_object = QtGui.QApplication(sys.argv)

	simgui = SimGui()
	simgui.show()
	
	application_object.exec_()


if __name__=="__main__":
        main(sys.argv)
