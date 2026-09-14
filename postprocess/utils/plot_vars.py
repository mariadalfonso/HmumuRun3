"""
plot_vars.py -- lookup functions for plotting variables.

Contains only the variables histo_config.get_active_vars() actually plots.
See config/branches.yaml for what the snapshots make available.
"""


def get_expr(varname):
    """C++ expression for this variable, evaluated by RDataFrame.Define()."""
    exprs = {
        # --- mass ---
        'dimu_mass': 'HiggsCandCorrMass',
        # --- discriminants ---
        'mva': 'discrMVA0',
        'category_vlcat': 'category',
        'category_ttlcat': 'category',
        'category_tthcat': 'category',
        # --- muons / dimuon ---
        'muon1_pt': 'Muon1_pt',
        'muon2_pt': 'Muon2_pt',
        'muon1_eta': 'Muon1_eta',
        'muon2_eta': 'Muon2_eta',
        'muon1_norm_pt': 'Muon1_pt/HiggsCandCorrMass',
        'muon2_norm_pt': 'Muon2_pt/HiggsCandCorrMass',
        'muon1_sip3d': 'Muon1_sip3d',
        'muon2_sip3d': 'Muon2_sip3d',
        'deta_muons': 'fabs(Muon1_eta-Muon2_eta)',
        'dimuon_pt': 'HiggsCandCorrPt',
        'dimuon_eta': 'HiggsCandCorrEta',
        'dimuon_rapidity': 'HiggsCandCorrRapidity',
        'costhetacs': 'cosThetaCS',
        'phistarcs': 'phiStarCS',
        # --- jets / MET ---
        'njets': 'nGoodJetsAll',
        'jet1_pt': 'Jet1_Pt',
        'jet1_eta': 'Jet1_Eta',
        'met_pt': 'PuppiMET_pt',
    }
    return exprs.get(varname, varname)  # fallback: assume the name is the branch


def get_binning(varname):
    """(nbins, low, high). None means no binning is defined for this variable,
    and makeHistos.py skips it."""
    binning = {
        # --- mass ---
        'dimu_mass': (130, 70.0, 200.0),
        # --- discriminants ---
        'mva': (100, 0.0, 1.0),
        'category_vlcat': (5, 0.0, 5.0),
        'category_ttlcat': (6, 0.0, 6.0),
        'category_tthcat': (4, 0.0, 4.0),
        # --- muons / dimuon ---
        'muon1_pt': (100, 0.0, 200.0),
        'muon2_pt': (100, 0.0, 200.0),
        'muon1_eta': (60, -3.0, 3.0),
        'muon2_eta': (60, -3.0, 3.0),
        'muon1_norm_pt': (100, 0.0, 2.0),
        'muon2_norm_pt': (100, 0.0, 2.0),
        'muon1_sip3d': (100, 0.0, 20.0),
        'muon2_sip3d': (100, 0.0, 20.0),
        'deta_muons': (60, 0.0, 6.0),
        'dimuon_pt': (100, 0.0, 400.0),
        'dimuon_eta': (60, -3.0, 3.0),
        'dimuon_rapidity': (60, -3.0, 3.0),
        'costhetacs': (50, -1.2, 1.2),
        'phistarcs': (50, -4, 4),
        # --- jets / MET ---
        'njets': (10, 0.0, 10.0),
        'jet1_pt': (100, 0.0, 300.0),
        'jet1_eta': (100, -5.0, 5.0),
        'met_pt': (100, 0.0, 300.0),
    }
    return binning.get(varname)


def get_xlabel(varname):
    """x-axis title (ROOT TLatex). Falls back to the variable name."""
    labels = {
        # --- mass ---
        'dimu_mass': 'm_{#mu#mu} [GeV]',
        # --- discriminants ---
        'mva': 'MVA discr',
        'category_vlcat': 'category',
        'category_ttlcat': 'category',
        'category_tthcat': 'category',
        # --- muons / dimuon ---
        'muon1_pt': 'p_{T}^{#mu_{1}} [GeV]',
        'muon2_pt': 'p_{T}^{#mu_{2}} [GeV]',
        'muon1_eta': '#eta^{#mu_{1}}',
        'muon2_eta': '#eta^{#mu_{2}}',
        'muon1_norm_pt': 'p_{T}^{#mu_{1}}/m_{#mu#mu}',
        'muon2_norm_pt': 'p_{T}^{#mu_{2}}/m_{#mu#mu}',
        'muon1_sip3d': 'Muon1_sip3d',
        'muon2_sip3d': 'Muon2_sip3d',
        'deta_muons': '|#Delta#eta(#mu_{1}, #mu_{2})|',
        'dimuon_pt': 'p_{T}^{#mu#mu} [GeV]',
        'dimuon_eta': '#eta_{#mu#mu}',
        'dimuon_rapidity': 'y_{#mu#mu}',
        'costhetacs': 'cos#theta^{*}_{CS}',
        'phistarcs': '#phi^{*}_{CS}',
        # --- jets / MET ---
        'njets': 'Njets',
        'jet1_pt': 'Jet1 p_{T} [GeV]',
        'jet1_eta': 'Jet1 #eta',
        'met_pt': 'PuppiMET p_{T} [GeV]',
    }
    return labels.get(varname, varname)


def get_logy_vars():
    """Variables plotted with a log y-axis. Opt-in list; anything not here
    is drawn linear."""
    return [
        'dimu_mass',
        'mva',
        'muon1_pt',
        'muon2_pt',
        'muon1_eta',
        'muon2_eta',
        'dimuon_pt',
        'muon1_sip3d',
        'muon2_sip3d',
        'deta_muons',
        'njets',
    ]


def get_bin_labels(varname):
    """Custom x-axis bin labels for categorical (category-index) variables."""
    bin_labels = {
        "category_vlcat": ["", "H_{#mu#mu}+e", "H_{#mu#mu}+ee", "H_{#mu#mu}+#mu",
                            "H_{#mu#mu}+#mu#mu", "H_{#mu#mu}+e#mu"],
        "category_ttlcat": ["", "H_{#mu#mu}+e", "H_{#mu#mu}+ee", "H_{#mu#mu}+#mu",
                             "H_{#mu#mu}+#mu#mu", "H_{#mu#mu}+e#mu"],
        "category_tthcat": ["", "n_top", "n_W", "resolved (5jets)"],
    }
    return bin_labels.get(varname)


def get_ratio_range(varname):
    """(low, high) for the ratio-pad y-axis. Falls back to a wide default."""
    ranges = {
        "dimu_mass": (0.90, 1.10),
    }
    return ranges.get(varname, (0.5, 1.5))


def get_blind_range(varname):
    """(low, high) mass window blinded in data, if this variable is blinded."""
    blinds = {
        "dimu_mass": (110, 150),
    }
    return blinds.get(varname)
