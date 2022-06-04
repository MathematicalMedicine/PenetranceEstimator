"""GUI front-end to the Penetrance Estimator app."""

from pathlib import Path
from ast import literal_eval
from enum import Enum

from PyQt5 import uic, QtWidgets as QtW, QtCore
import matplotlib
import seaborn as sns
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

import penEst as pe
from sinaplot import sinaplot

#matplotlib.rcParams["text.usetex"] = True
        # turns out this was unnecessary, woo

# Canvas and MPL widget management inspired by/borrowed from
# https://stackoverflow.com/a/44029435

class PenEstCanvas(FigureCanvasQTAgg):
    """matplotlib canvas and figure for the Penetrance Estimator app."""
    
    def __init__(self):
        self.figure = Figure()
        self.figure.set_tight_layout(True)
        self.ax_ftilde = self.figure.add_subplot(1, 2, 1)
        self.ax_ftildestar = self.figure.add_subplot(1, 2, 2)
        FigureCanvasQTAgg.__init__(self, self.figure)
        self.setSizePolicy(QtW.QSizePolicy.Expanding, QtW.QSizePolicy.Expanding)
        self.updateGeometry()
    
    def plot_data(self, df):
        """Given our dataframe, clears our axes and plots the new data. """
        
        self.ax_ftilde.clear()
        self.ax_ftildestar.clear()
        
        # Are we Figure 1 or Figure 2?
        try:
            df['SimNum']
        except KeyError:
            # Lineplot (figure 2)
            sns.lineplot(data=df.pivot(pe.ALPHA, "k", pe.FT),
                    ax=self.ax_ftilde, markers=True)
            sns.lineplot(data=df.pivot(pe.ALPHA, "k", pe.FTS),
                    ax=self.ax_ftildestar, markers=True)
        else:
            # Sinaplot (figure 1)
            sinaplot(x="N", y=pe.FT, data=df, ax=self.ax_ftilde,
                    edgecolor="black", alpha=.5, violin=False)
            sinaplot(x="N", y=pe.FTS, data=df, ax=self.ax_ftildestar,
                    edgecolor="black", alpha=.5, violin=False)
        
        self.ax_ftilde.set_ylabel(pe.FT, rotation=0)
        self.ax_ftildestar.set_ylabel(pe.FTS, rotation=0)
        self.ax_ftilde.set_ylim(ymin=0, ymax=1)
        self.ax_ftildestar.set_ylim(ymin=0, ymax=1)
        self.ax_ftilde.axhline(0.5, alpha=0.5, dashes=(5,2))
        self.ax_ftildestar.axhline(0.5, alpha=0.5, dashes=(5,2))
        
        self.draw()

class PenEstPlot(QtW.QWidget):
    """Widget representing the space for a matplotlib canvas."""
    
    def __init__(self, parent):
        super().__init__(parent)
        # FIXME: at this point, remove the stylesheet
        self.canvas = PenEstCanvas()
        self.vbl = QtW.QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.setLayout(self.vbl)


class PenEstFigureGroup(QtW.QGroupBox):
    """Widget representing the grouping of a plot and all its input widgets."""
    
    def __init__(self, parent):
        super().__init__(parent)
        self.dataframe = None
    
    def widgetcheck(self):
        """Pulls all the widgets we need from the main window for this plot."""
        # This is a separate method because there's some sort of strange race
        # condition involved with self.objectName() such that it doesn't seem
        # to be set while __init__ is running. So instead of doing this at init
        # time, I do this the first time we refresh the plot.
        
        try:
            self.PlotWidget
        except AttributeError:
            # Pulling our relevant widgets, since they're not automatically
            # added as children
            widgetstore = self.parent().parent()
            name = self.objectName()
            
            for widgettype in ("PlotWidget", "VarAlpha", "VarBeta", "VarS",
                    "VarK", "VarRep"):
                try:
                    setattr(self, widgettype,
                            getattr(widgetstore, f"{name}_{widgettype}"))
                except AttributeError:
                    pass
    
    def newdataframe(self):
        """Returns a new dataframe built from the user-provided inputs."""
        self.widgetcheck()
        
        # alpha - If user-specifiable, it's one value. Otherwise, it's a fixed
        # list of values.
        try:
            alpha = self.VarAlpha.value()
        except AttributeError:
            alpha = [0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5]
        # beta = Always user-specified with the same type of control
        beta = self.VarBeta.value()
        # s - Also always user-specified with the same type of control
        s = self.VarS.value()
        # k - Either a single value, OR a list of values
        try:
            k = self.VarK.value()
        except AttributeError:
            k = literal_eval(self.VarK.text())
        # N - A fixed list of values only if NumSims exists.
        # NumSims - Either it's user-specified or it doesn't exist.
        try:
            NumSims = self.VarRep.value()
        except AttributeError:
            NumSims = None
            N = None
        else:
            N = [1, 2, 4, 6, 8, 10, 30, 50]
        
        progress = ProxiedProgressDialog(self.parent().parent(),
                labeltext="Calculating figure data. Please wait...")
        return pe.calc_stats(alpha, beta, s, k, N, NumSims, progress=progress)
    
    def refreshPlot(self):
        """Plots stats for the current entered variable values."""
        
        self.dataframe = self.newdataframe()
        self.PlotWidget.canvas.figure.suptitle(
                self.title()[10:].replace("𝛂", pe.ALPHA))
        self.PlotWidget.canvas.plot_data(self.dataframe)
    
    def savePlotImage(self):
        """Saves our plot to an image file."""
        
        filterstring = ' '.join([f"*.{ext}" for ext in
                self.PlotWidget.canvas.get_supported_filetypes().keys()])
        fileselect = QtW.QFileDialog.getSaveFileName(self,
                "Save Image", filter=f"Images ({filterstring})")
        
        if fileselect[0]:  #verify cancel wasn't pressed
            self.PlotWidget.canvas.figure.savefig(fileselect[0])
    
    def saveStatsFile(self):
        """Saves the current plotted data to a text file."""
        
        fileselect = QtW.QFileDialog.getSaveFileName(self,
                "Save Stats", filter="Text files (*.txt)")
        
        if fileselect[0]:
            with open(fileselect[0], 'w') as outfile:
                with pe.DATA_PRINT_OPTIONS:
                    outfile.write(self.dataframe)


UI_MainWindow = uic.loadUiType(str(
        Path(__file__).parent.joinpath("penestapp.ui")))[0]
class PenEstAppMain(QtW.QMainWindow, UI_MainWindow):
    """Main (only) window for a GUI program for calculating and plotting
    our penetrance estimations.
    
    """
    
    def __init__(self):
        super().__init__()
        self.setupUi(self)
        self.setWindowTitle("Penetrance Estimator")
        
        # showing an initial plot (using default values that are embedded in
        # penestapp.ui) so there's something shinyesque to look at
        self.Fig1.refreshPlot()
        self.Fig2.refreshPlot()
    
    # Should I possibly add a "Save Combined Plot" button here?
    # Would be basically pe.generate_and_display_test_figure with plt.savefig()
    # instead of plt.show()...


class ProxiedProgressDialog(QtCore.QObject):
    """A very simple proxy object setup for a QProgressDialog that sends it
    updates via signals and slots rather than direct method calls (because
    apparently trying to do direct method calls in a worker thread is unsafe).
    
    """
    
    # Signals needed for communication with our progress dialog.
    showtime = QtCore.Signal()
    incr_progress = QtCore.Signal(int)
    set_max_progress = QtCore.Signal(int)
    task_complete = QtCore.Signal()
    
    def __init__(self, mainwindow, labeltext="[progress]"):
        super().__init__()
        
        # Track the state of our progress. Needed to determine when we actually
        # close the progress dialog.
        self.progress_state = 0
        self.max_progress = 1
        self._labeltext = labeltext
        
        # Set up progress dialog.
        # We explicitly do NOT save a reference to same because that would
        # invite temptation to call methods on it directly, which would defeat
        # the whole purpose of having this here proxy object - can't manipulate
        # the dialog directly within a worker thread. Instead, we make sure we
        # have signals for every method we might want to call, and use those
        # methods as the associated slots.
        self.progdlg = QtW.QProgressDialog(mainwindow)
        self.progdlg.setWindowTitle("Calculating...")
        self.message = labeltext
        self.progdlg.setCancelButton(None) # FIXME: make a working cancel button
        self.progdlg.setWindowModality(QtCore.Qt.WindowModal)
        
        # Slot connections for all QProgressDialog methods we care about.
        self.showtime.connect(self.progdlg.show)
        self.incr_progress.connect(self.progdlg.setValue)
        self.set_max_progress.connect(self.progdlg.setMaximum)
        self.task_complete.connect(self.progdlg.close)
    
    def tick(self):
        """Increments progress state, and emits a signal to the progress dialog
        to do the same. Automatically signals the progress dialog to close when
        we increment beyond our maximum progress.
        
        """
        
        self.progress_state += 1
        if self.progress_state <= self.max_progress:
            self.incr_progress.emit(self.progress_state)
        else:
            self.task_complete.emit()
    
    @property
    def total(self):
        """Returns our max progress value."""
        return self.max_progress
    
    @total.setter
    def total(self, max_progress):
        """Informs the progress dialog as to what our actual max value is, then
        starts it.
        
        """
        
        self.max_progress = max_progress
        self.set_max_progress.emit(max_progress)
        self.showtime.emit()
    
    @property
    def message(self):
        """Returns our current label text."""
        return self._labeltext
    
    @message.setter
    def message(self, message):
        """Updates our progress bar status message."""
        self._labeltext = message
        # Edges of the label get cut off on later versions of macOS, because
        # reasons, so adding padding.
        self.progdlg.setLabelText(f"    {message}    ")
        
    
    def killme(self):
        self.task_complete.emit()
