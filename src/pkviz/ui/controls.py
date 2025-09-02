# src/pkviz/ui/controls.py
from dataclasses import dataclass, field
from PySide6.QtCore import Signal, QObject
from PySide6.QtWidgets import QWidget, QVBoxLayout, QPushButton, QDoubleSpinBox, QSpinBox, QComboBox, QFrame, QLabel
import math
from pkengine.types import Drug, Regimen
from pkengine.dosing import fixed_every_n_days, single_dose, iv_infusion

@dataclass
class SimulateRequest:
    drugs: list[Drug] = field(default_factory=list)
    regimen: Regimen | None = None
    t_end_h: float = 8*7*24
    dt_h: float = 1.0

class ControlsPanel(QFrame):
    simulateRequested = Signal(SimulateRequest)

    def __init__(self):
        super().__init__()
        self.setFrameShape(QFrame.StyledPanel)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Controls"))

        # --- Dosing Parameters ---
        layout.addWidget(QLabel("Dosing"))
        self.dose = QDoubleSpinBox(); self.dose.setRange(1, 1e6); self.dose.setValue(250)
        self.dose.setSuffix(" mg")
        layout.addWidget(QLabel("Dose (mg)"))
        layout.addWidget(self.dose)
        
        self.interval = QSpinBox(); self.interval.setRange(1, 60); self.interval.setValue(3)
        self.interval.setSuffix(" days")
        layout.addWidget(QLabel("Interval (days)"))
        layout.addWidget(self.interval)
        
        self.weeks = QSpinBox(); self.weeks.setRange(1, 104); self.weeks.setValue(8)
        self.weeks.setSuffix(" weeks")
        layout.addWidget(QLabel("Duration (weeks)"))
        layout.addWidget(self.weeks)
        
        self.route = QComboBox(); self.route.addItems(["im", "iv_bolus", "iv_infusion"])
        layout.addWidget(QLabel("Route"))
        layout.addWidget(self.route)

        # Infusion-only parameter: duration (hours)
        self.infusion_duration = QDoubleSpinBox(); self.infusion_duration.setDecimals(2)
        self.infusion_duration.setRange(0.1, 1e5); self.infusion_duration.setValue(2.0)
        self.infusion_duration.setSuffix(" h")
        self.lbl_infusion_duration = QLabel("Infusion duration (h)")
        layout.addWidget(self.lbl_infusion_duration)
        layout.addWidget(self.infusion_duration)
        # Show duration field only when route is iv_infusion
        self.route.currentIndexChanged.connect(self._update_route_visibility)
        self._update_route_visibility()
        
        # --- PK Parameters ---
        layout.addWidget(QLabel("PK Parameters"))

        # Mode selector
        layout.addWidget(QLabel("Input mode"))
        self.mode = QComboBox(); self.mode.addItems(["Clinical (t½, V, t½a)", "Mechanistic (CL, V, ka)"])
        layout.addWidget(self.mode)

        # Shared parameter: V
        self.V = QDoubleSpinBox(); self.V.setRange(1, 20000); self.V.setValue(300.0)
        self.V.setSuffix(" L")
        self.lbl_V = QLabel("V (volume, L)")
        layout.addWidget(self.lbl_V)
        layout.addWidget(self.V)

        # Mechanistic inputs: CL, ka
        self.ka = QDoubleSpinBox(); self.ka.setDecimals(3); self.ka.setRange(0.001, 10); self.ka.setValue(0.10)
        self.ka.setSuffix(" /h")
        self.lbl_ka = QLabel("ka (absorption rate, /h)")
        layout.addWidget(self.lbl_ka)
        layout.addWidget(self.ka)

        self.CL = QDoubleSpinBox(); self.CL.setRange(0.1, 5000); self.CL.setValue(12.0)
        self.CL.setSuffix(" L/h")
        self.lbl_CL = QLabel("CL (clearance, L/h)")
        layout.addWidget(self.lbl_CL)
        layout.addWidget(self.CL)

        # Clinical inputs: t_half (elim) and t_half_abs (ka)
        self.t_half = QDoubleSpinBox(); self.t_half.setDecimals(2); self.t_half.setRange(0.1, 1e5)
        self.t_half.setSuffix(" h")
        self.lbl_t_half = QLabel("t½ (elimination half-life, h)")
        layout.addWidget(self.lbl_t_half)
        layout.addWidget(self.t_half)

        self.t_half_abs = QDoubleSpinBox(); self.t_half_abs.setDecimals(2); self.t_half_abs.setRange(0.1, 1e5)
        self.t_half_abs.setSuffix(" h")
        self.lbl_t_half_abs = QLabel("t½a (absorption half-life, h)")
        layout.addWidget(self.lbl_t_half_abs)
        layout.addWidget(self.t_half_abs)

        # Initialize clinical values from mechanistic defaults
        ln2 = math.log(2.0)
        self.t_half.setValue((ln2 * float(self.V.value())) / float(self.CL.value()))
        self.t_half_abs.setValue(ln2 / float(self.ka.value()))

        # Default to Clinical mode, set visibility
        self.mode.setCurrentIndex(0)
        self.mode.currentIndexChanged.connect(self._update_mode_visibility)
        self._update_mode_visibility()

        # Sampling parameter: dt (hours)
        self.dt = QDoubleSpinBox(); self.dt.setDecimals(2)
        self.dt.setRange(0.01, 24.0); self.dt.setValue(1.0)
        self.dt.setSuffix(" h")
        self.lbl_dt = QLabel("Sampling dt (h)")
        layout.addWidget(self.lbl_dt)
        layout.addWidget(self.dt)

        go = QPushButton("Simulate"); layout.addWidget(go)
        go.clicked.connect(self._emit_request)

    def _update_mode_visibility(self):
        clinical = (self.mode.currentIndex() == 0)
        # Clinical widgets visible?
        self.lbl_t_half.setVisible(clinical)
        self.t_half.setVisible(clinical)
        self.lbl_t_half_abs.setVisible(clinical)
        self.t_half_abs.setVisible(clinical)
        # Mechanistic widgets visible?
        self.lbl_CL.setVisible(not clinical)
        self.CL.setVisible(not clinical)
        self.lbl_ka.setVisible(not clinical)
        self.ka.setVisible(not clinical)

    def _update_route_visibility(self):
        is_infusion = (self.route.currentText() == "iv_infusion")
        self.lbl_infusion_duration.setVisible(is_infusion)
        self.infusion_duration.setVisible(is_infusion)

    def _emit_request(self):
        # build Drug
        ln2 = math.log(2.0)
        if self.mode.currentIndex() == 0:  # Clinical mode
            V = float(self.V.value())
            t_half = float(self.t_half.value())
            t_half_abs = float(self.t_half_abs.value())
            CL = (ln2 * V) / t_half
            ka = ln2 / t_half_abs
        else:  # Mechanistic mode
            V = float(self.V.value())
            CL = float(self.CL.value())
            ka = float(self.ka.value())

        drug = Drug(
            drug_id="drug",  # later make editable
            CL_L_per_h=CL,
            V_L=V,
            ka_per_h=ka,
        )
        # build Regimen (start with simple qNd IM; add route branching next)
        route = self.route.currentText()
        if route == "iv_bolus":
            reg = single_dose(amount_mg=float(self.dose.value()), start_h=0.0, route="iv_bolus", drug_id=drug.drug_id)
        elif route == "iv_infusion":
            reg = iv_infusion(amount_mg=float(self.dose.value()), start_h=0.0, duration_h=float(self.infusion_duration.value()), drug_id=drug.drug_id)
        else:
            reg = fixed_every_n_days(amount_mg=float(self.dose.value()),
                                     every_days=int(self.interval.value()),
                                     weeks=int(self.weeks.value()),
                                     route="im", drug_id=drug.drug_id)
        t_end_h = int(self.weeks.value()) * 7 * 24
        req = SimulateRequest(drugs=[drug], regimen=reg, t_end_h=t_end_h, dt_h=float(self.dt.value()))
        self.simulateRequested.emit(req)