from collections import defaultdict
from .types import Regimen, Dose

def split_regimen_by_drug(regimen: Regimen) -> dict[str, Regimen]:
    """
    Group doses by drug_id into separate Regimens.
    """
    buckets: dict[str, list[Dose]] = defaultdict(list)
    for d in regimen.doses:
        buckets[d.drug_id].append(d)
    return {
        drug_id: Regimen(doses=tuple(sorted(ds, key=lambda x: (x.start_h, x.route))))
        for drug_id, ds in buckets.items()
    }