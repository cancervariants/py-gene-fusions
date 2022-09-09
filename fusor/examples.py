"""Provide programmatic access to example objects."""
import json

from fusor import APP_ROOT
from fusor.models import AssayedFusion, CategoricalFusion

examples_dir = APP_ROOT / "examples"

with open(examples_dir / "alk.json") as f:
    alk = CategoricalFusion(**json.load(f))

with open(examples_dir / "bcr_abl1.json") as f:
    bcr_abl1 = CategoricalFusion(**json.load(f))

with open(examples_dir / "bcr_abl1_expanded.json") as f:
    bcr_abl1_expanded = CategoricalFusion(**json.load(f))

with open(examples_dir / "ewsr1.json") as f:
    ewsr1 = AssayedFusion(**json.load(f))

with open(examples_dir / "igh_myc.json") as f:
    igh_myc = CategoricalFusion(**json.load(f))

with open(examples_dir / "tpm3_ntrk1.json") as f:
    tpm3_ntrk1 = AssayedFusion(**json.load(f))

with open(examples_dir / "tpm3_pdgfrb.json") as f:
    tpm3_pdgfrb = AssayedFusion(**json.load(f))
