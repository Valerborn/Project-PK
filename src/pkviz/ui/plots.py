# src/pkviz/ui/plots.py
from PySide6.QtWidgets import QWidget, QVBoxLayout
import pyqtgraph as pg


class PlotWidget(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)

        # Main plot area
        self.plot_widget = pg.PlotWidget()
        self.plot_widget.setLabel("left", "Concentration", units="mg/L")
        self.plot_widget.setLabel("bottom", "Time", units="h")
        self.plot_widget.showGrid(x=True, y=True, alpha=0.3)
        self.plot_widget.addLegend()
        layout.addWidget(self.plot_widget)

        self.curves = {}  # store references for updates

    def plot_curves(self, results: dict[str, tuple]):
        self.plot_widget.clear()
        self.curves = {}
        for label, (t, C) in results.items():
            curve = self.plot_widget.plot(
                t, C,
                pen=pg.mkPen(width=2),
                name=label
            )
            self.curves[label] = curve

    def clear(self):
        self.plot_widget.clear()
        self.curves = {}