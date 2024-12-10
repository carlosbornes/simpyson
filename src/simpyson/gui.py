# src/simpyson/gui.py

from PyQt5 import QtWidgets
import sys
from PyQt5.QtWidgets import (
    QFileDialog, QMessageBox, QVBoxLayout, QWidget, QMainWindow, QAction, QInputDialog,
    QSplitter, QListWidget, QHBoxLayout
)
from PyQt5.QtCore import Qt
from PyQt5.QtWebEngineWidgets import QWebEngineView
import plotly.graph_objects as go
from simpyson.io import SimpReader
import numpy as np
import os

class SimpysonGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Simpyson GUI")
        self.setGeometry(100, 100, 1200, 800)

        self.files_data = {}  # Dictionary to store {filename: data}
        self.current_file = None
        self.data = None
        self.filename = None  # Added to store full path

        self.init_ui()

    def init_ui(self):
        self.create_menu()
        self.create_main_layout()

    def create_main_layout(self):
        # Create central widget with horizontal layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QHBoxLayout(central_widget)

        # Create splitter for resizable panels
        splitter = QSplitter(Qt.Horizontal)
        
        # Create file list widget with multiple selection
        self.file_list = QListWidget()
        self.file_list.setSelectionMode(QListWidget.ExtendedSelection)
        self.file_list.itemSelectionChanged.connect(self.on_selection_changed)
        splitter.addWidget(self.file_list)

        # Create plot area
        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)
        self.browser = QWebEngineView()
        plot_layout.addWidget(self.browser)
        plot_widget.setLayout(plot_layout)
        
        splitter.addWidget(plot_widget)
        splitter.setSizes([200, 1000])  # Left panel 200px, rest for plot
        
        layout.addWidget(splitter)

    def create_menu(self):
        menubar = self.menuBar()

        # File Menu
        file_menu = menubar.addMenu('File')

        open_action = QAction('Open', self)
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)

        save_action = QAction('Save', self)
        save_action.triggered.connect(self.save_file)
        file_menu.addAction(save_action)

        file_menu.addSeparator()

        exit_action = QAction('Exit', self)
        exit_action.triggered.connect(QtWidgets.qApp.quit)
        file_menu.addAction(exit_action)

        # Process Menu  
        process_menu = menubar.addMenu('Process')

        convert_hz_to_ppm = QAction('Convert Hz to ppm', self)
        convert_hz_to_ppm.triggered.connect(self.convert_hz_to_ppm)
        process_menu.addAction(convert_hz_to_ppm)

        convert_fid_to_spe = QAction('Convert FID to SPE', self)
        convert_fid_to_spe.triggered.connect(self.convert_fid_to_spe)
        process_menu.addAction(convert_fid_to_spe)

        convert_spe_to_fid = QAction('Convert SPE to FID', self)
        convert_spe_to_fid.triggered.connect(self.convert_spe_to_fid)
        process_menu.addAction(convert_spe_to_fid)

    def save_file(self):
        if not self.current_file:
            QMessageBox.warning(self, 'Save File', 'No file selected!')
            return

        if self.data:
            options = QFileDialog.Options()
            save_filename, _ = QFileDialog.getSaveFileName(
                self, 'Save File', '', 'CSV Files (*.csv)', options=options
            )
            if save_filename:
                if 'hz' in self.data.data:
                    x_data = self.data.data['hz']
                    x_label = 'Hz'
                elif 'ppm' in self.data.data:
                    x_data = self.data.data['ppm']
                    x_label = 'ppm'
                elif 'time' in self.data.data:
                    x_data = self.data.data['time']
                    x_label = 'Time'
                else:
                    QMessageBox.warning(self, 'Save File', 'Data does not contain hz, ppm, or time.')
                    return

                try:
                    np.savetxt(
                        save_filename,
                        np.column_stack((x_data, self.data.data['real'])),
                        delimiter=",",
                        header=f"{x_label},Real",
                        comments=""
                    )
                    QMessageBox.information(self, 'Save File', 'File saved successfully!')
                except Exception as e:
                    QMessageBox.warning(self, 'Save File', f'Error saving file: {str(e)}')
            else:
                QMessageBox.warning(self, 'Save File', 'No filename provided!')
        else:
            QMessageBox.warning(self, 'Save File', 'No data to save!')

    def open_file(self):
        options = QFileDialog.Options()
        filenames, _ = QFileDialog.getOpenFileNames(
            self, 'Open File', '', 'SIMPSON Files (*.spe *.fid)', options=options
        )
        
        for filename in filenames:
            if filename:
                file_format = filename.split('.')[-1]
                base_name = os.path.basename(filename)
                
                # Store both the data and full path
                self.files_data[base_name] = {
                    'data': SimpReader(filename, format=file_format),
                    'path': filename
                }
                
                # Add to list widget if not already there
                if self.file_list.findItems(base_name, Qt.MatchExactly) == []:
                    self.file_list.addItem(base_name)
                
                # If this is the first file, display it
                if not self.current_file:
                    self.current_file = base_name
                    self.data = self.files_data[base_name]['data']
                    self.filename = filename
                    self.plot_data()

    def on_selection_changed(self):
        selected_items = self.file_list.selectedItems()
        if selected_items:
            self.current_file = selected_items[0].text()
            self.data = self.files_data[self.current_file]['data']
            self.filename = self.files_data[self.current_file]['path']
            self.plot_data(selected_items)

    def plot_data(self, selected_items=None):
        if not selected_items:
            selected_items = [self.file_list.findItems(self.current_file, Qt.MatchExactly)[0]] if self.current_file else []
        
        if selected_items:
            fig = go.Figure()
            
            # Determine x-axis type from first selected item
            first_data = self.files_data[selected_items[0].text()]['data']
            if 'ppm' in first_data.data:
                x_axis = 'ppm'
                xlabel = 'Chemical Shift (ppm)'
            elif 'hz' in first_data.data:
                x_axis = 'hz'
                xlabel = 'Frequency (Hz)'
            elif 'time' in first_data.data:
                x_axis = 'time'
                xlabel = 'Time (ms)'
            else:
                QMessageBox.warning(self, 'Plot Data', 'Data does not contain valid axis information.')
                return
    
            # Plot each selected spectrum
            for item in selected_items:
                data = self.files_data[item.text()]['data']
                if x_axis in data.data:
                    fig.add_trace(go.Scatter(
                        x=data.data[x_axis],
                        y=data.data['real'],
                        name=item.text(),
                        mode='lines'
                    ))
    
            # Update layout with consistent axis settings
            fig.update_layout(
                title='NMR Spectrum',
                xaxis_title=xlabel,
                yaxis_title='Intensity',
                showlegend=True
            )
    
            # Always invert x-axis for ppm and hz
            if x_axis in ['ppm', 'hz']:
                fig.update_xaxes(autorange="reversed")
    
            # Update plot display
            html_content = fig.to_html(include_plotlyjs='cdn', full_html=True)
            self.browser.setHtml(html_content)
        else:
            QMessageBox.warning(self, 'Plot Data', 'No data to plot!')

    def convert_hz_to_ppm(self):
        if not self.current_file:
            QMessageBox.warning(self, 'Convert Hz to ppm', 'No file selected!')
            return
            
        if self.data and 'hz' in self.data.data:
            b0, ok = QInputDialog.getText(self, 'Input', 'Enter B0 (e.g., 400MHz or 9.4T):')
            if ok and b0:
                nucleus, ok = QInputDialog.getText(self, 'Input', 'Enter nucleus (e.g., 1H or 13C):')
                if ok and nucleus:
                    try:
                        # Re-read the file with b0 and nucleus values
                        file_format = self.filename.split('.')[-1]
                        new_data = SimpReader(self.filename, format=file_format, b0=b0, nucleus=nucleus)
                        self.files_data[self.current_file]['data'] = new_data
                        self.data = new_data
                        self.plot_data()
                    except ValueError as e:
                        QMessageBox.warning(self, 'Convert Hz to ppm', str(e))
                else:
                    QMessageBox.warning(self, 'Convert Hz to ppm', 'Nucleus must be specified!')
            else:
                QMessageBox.warning(self, 'Convert Hz to ppm', 'B0 must be specified!')
        else:
            QMessageBox.warning(self, 'Convert Hz to ppm', 'No data to convert!')

    def convert_fid_to_spe(self):
        if not self.current_file:
            QMessageBox.warning(self, 'Convert FID to SPE', 'No file selected!')
            return
            
        if self.data and self.data.format == 'fid':
            try:
                new_data = self.data.to_spe()
                self.files_data[self.current_file]['data'] = new_data
                self.data = new_data
                self.plot_data()
            except ValueError as e:
                QMessageBox.warning(self, 'Convert FID to SPE', str(e))
        else:
            QMessageBox.warning(self, 'Convert FID to SPE', 'No FID data to convert!')

    def convert_spe_to_fid(self):
        if not self.current_file:
            QMessageBox.warning(self, 'Convert SPE to FID', 'No file selected!')
            return
            
        if self.data and self.data.format == 'spe':
            try:
                new_data = self.data.to_fid()
                self.files_data[self.current_file]['data'] = new_data
                self.data = new_data
                self.plot_data()
            except ValueError as e:
                QMessageBox.warning(self, 'Convert SPE to FID', str(e))
        else:
            QMessageBox.warning(self, 'Convert SPE to FID', 'No SPE data to convert!')

def main():
    app = QtWidgets.QApplication(sys.argv)
    gui = SimpysonGUI()
    gui.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()