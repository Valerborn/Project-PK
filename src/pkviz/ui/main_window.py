# src/pkviz/ui/main_window.py (skeleton essentials)
from PySide6.QtWidgets import QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QStatusBar
from .controls import ControlsPanel, SimulateRequest
from .plots import PlotWidget
from pkengine.simulate import run_single, run_multi
from pkengine.metrics import cmax, tmax, auc_trapz, peak_to_trough_ratio, fluctuation_index

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("PK Viz")
        self.resize(1100, 680)

        central = QWidget(self); self.setCentralWidget(central)
        root = QHBoxLayout(central)

        self.controls = ControlsPanel()
        self.plot = PlotWidget()
        root.addWidget(self.controls, 0)
        root.addWidget(self.plot, 1)

        self.status = QStatusBar(); self.setStatusBar(self.status)

        # wire events
        self.controls.simulateRequested.connect(self.on_simulate)

        # first run using current control values
        self.controls._emit_request()

    def on_simulate(self, req: SimulateRequest):
        try:
            # single vs multi (based on req.drugs list)
            if len(req.drugs) == 1:
                drug = req.drugs[0]
                t, C = run_single(drug, req.regimen, t_end_h=req.t_end_h, dt_h=req.dt_h)
                self.plot.plot_curves({drug.drug_id: (t, C)})
                # quick metrics example:
                msg = f"Cmax {cmax(C):.3f} mg/L at {tmax(t, C):.1f} h | AUC {auc_trapz(t, C):.1f}"
                self.status.showMessage(msg, 5000)
            else:
                results = run_multi({d.drug_id: d for d in req.drugs}, req.regimen,
                                    t_end_h=req.t_end_h, dt_h=req.dt_h)
                self.plot.plot_curves(results)
                self.status.showMessage(f"Simulated {len(results)} drugs", 5000)
        except Exception as e:
            self.status.showMessage(f"Error: {e}", 8000)