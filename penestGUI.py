"""GUI front-end to the Penetrance Estimator app."""

from pathlib import Path
from ast import literal_eval
from enum import Enum

from PyQt5 import uic, QtWidgets as QtW
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg

import pe

class PlotType(Enum):
    numfam = 1
    alpha = 2

#matplotlib.rcParams["text.usetex"] = True
        # turns out this was unnecessary, woo

# Canvas and MPL widget management inspired by/borrowed from
# https://stackoverflow.com/a/44029435

class PenEstCanvas(FigureCanvasQTAgg):
    """matplotlib canvas and figure for the Penetrance Estimator app."""
    
    def __init__(self, plottype=PlotType.numfam):
        self.figure = Figure()
        self.ax_ftilde = self.figure.add_subplot(2, 1, 1)
        self.ax_ftildestar = self.figure.add_subplot(2, 1, 2)
        FigureCanvasQTAgg.__init__(self, self.figure)
        self.setSizePolicy(QtW.QSizePolicy.Expanding, QtW.QSizePolicy.Expanding)
        self.updateGeometry()
    
    def plot_data(self, df):
        """Given our dataframe, clears our axes and plots the new data. """
        
        self.ax_ftilde.clear()
        self.ax_ftildestar.clear()
        
        # Are we Figure 1 or Figure 2?
        if df['NumSims'][0] is not None:
            # Sinaplot (figure 1)
            sinaplot(x="N", y=pe.FT, data=df, ax=self.ax_ftilde,
                    edgecolor="black", alpha=.5, violin=False)
            sinaplot(x="N", y=pe.FTS, data=df, ax=self.ax_ftildestar,
                    edgecolor="black", alpha=.5, violin=False)
        else:
            # Lineplot (figure 2)
            sns.lineplot(data=df.pivot(pe.ALPHA, "k", pe.FT),
                    ax=self.ax_ftilde, markers=True)
            sns.lineplot(data=df.pivot(pe.ALPHA, "k", pe.FTS),
                    ax=self.ax_ftildestar, markers=True)
        
        self.ax_ftilde.set_ylabel(pe.FT, rotation=0)
        self.ax_ftildestar.set_ylabel(pe.FTS, rotation=0)
        self.ax_ftilde.set_ylim(ymin=0, ymax=1)
        self.ax_ftildestar.set_ylim(ymin=0, ymax=1)
        self.ax_ftilde.axhline(0.5, alpha=0.5, dashes=(5,2))
        self.ax_ftildestar.axhline(0.5, alpha=0.5, dashes=(5,2))
        
        self.tight_layout()
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
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dataframe = None
    
    @property
    def newdata(self):
        """Returns a new dataframe built from the user-provided inputs."""
        
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
        
        return pe.calc_stats(alpha, beta, s, k, N, NumSims)
    
    def refreshPlot(self):
        """Plots stats for the current entered variable values."""
        
        self.dataframe = self.newdata
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
        self.Figure1.refreshPlot()
        self.Figure2.refreshPlot()
    
    # Should I possibly add a "Save Combined Plot" button here?
    # Would be basically pe.generate_and_display_test_figure with plt.savefig()
    # instead of plt.show()...
