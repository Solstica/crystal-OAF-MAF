"""
Physical calibration target list.

Do not fill values by analogy with inorganic perovskites. Every numerical value
added here should be accompanied by a provenance entry in
`docs/model_notes/PARAMETER_PROVENANCE.md`.
"""

PARAMETER_TARGETS = {
    "crystal": [
        "local P-E / Landau-like constitutive coefficients",
        "background dielectric response",
        "gradient coefficient or domain-wall scale",
    ],
    "oaf": [
        "local susceptibility or P-E constitutive coefficients",
        "background dielectric response",
        "relaxation time or TDGL mobility proxy",
    ],
    "maf": [
        "dielectric susceptibility",
        "relaxation time or TDGL mobility proxy",
    ],
    "hhtt": [
        "mapping from HHTT content to phase fractions",
        "optional local switching-barrier correction if independently supported",
    ],
    "free_volume": [
        "PALS radius/fraction to dielectric response",
        "frequency dependence",
        "high-field constraint",
    ],
    "folded_boundary": [
        "boundary indicator to switching-barrier or strain coupling",
    ],
}
