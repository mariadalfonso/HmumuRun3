# pdfDefinitions.py
import ROOT
from ROOT import RooRealVar, RooDoubleCBFast, RooBernstein, RooChebychev
from ROOT import RooExponential, RooGaussian, RooAddPdf, RooGenericPdf, RooArgList, RooFormulaVar

class PDFDefinitions:
    """Centralized PDF definitions for signal and background fits"""
    
    @staticmethod
    def create_signal_pdfs(x, tag, suffix=""):
        """Create double Crystal Ball signal PDFs"""
        cb_mu = RooRealVar(f'cb_mu{suffix}{tag}', 'cb_mu', 125., 120., 130.)
        cb_sigma = RooRealVar(f'cb_sigma{suffix}{tag}', 'cb_sigma', 3, 0.5, 6.)
        cb_alphaL = RooRealVar(f'cb_alphaL{suffix}{tag}', 'cb_alphaL', 2., 0., 5.)
        cb_alphaR = RooRealVar(f'cb_alphaR{suffix}{tag}', 'cb_alphaR', 2., 0., 5.)
        cb_nL = RooRealVar(f'cb_nL{suffix}{tag}', 'cb_nL', 3., 0., 5.)
        cb_nR = RooRealVar(f'cb_nR{suffix}{tag}', 'cb_nR', 8., 0., 12.)

        pdf_cb = RooDoubleCBFast(f'crystal_ball{suffix}{tag}', 'crystal_ball', x, 
                                 cb_mu, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR)

        return {
            'pdf': pdf_cb,
            'params': {
                'mu': cb_mu, 'sigma': cb_sigma, 'alphaL': cb_alphaL,
                'alphaR': cb_alphaR, 'nL': cb_nL, 'nR': cb_nR
            }
        }

    @staticmethod
    def create_bkg_pdfs(x, tag):
        """
        Create all background PDFs - KEEP ALL OBJECTS IN MEMORY!
        
        Returns dict with:
            - 'pdfs': final PDF objects to use
            - 'params': all RooRealVar parameters
            - 'components': intermediate PDF components
            - 'formulas': all RooFormulaVar objects
        """
        
        params = {}
        pdfs = {}
        components = {}  # Store intermediate PDFs!
        formulas = {}    # Store formula variables!
        
        # ===== BERNSTEIN PDFs =====
        params['bern_c0'] = RooRealVar(f'bern_c0_{tag}', 'bern_c0', 1.8, 0.5, 2.)
        params['bern_c1'] = RooRealVar(f'bern_c1_{tag}', 'bern_c1', 1., 0., 2.)
        params['bern_c2'] = RooRealVar(f'bern_c2_{tag}', 'bern_c2', 0.5, 0., 1.5)
        params['bern_c3'] = RooRealVar(f'bern_c3_{tag}', 'bern_c3', 0.1, 0., 1.)
        params['bern_c4'] = RooRealVar(f'bern_c4_{tag}', 'bern_c4', 0.5, 0., 5.)
        params['bern_c5'] = RooRealVar(f'bern_c5_{tag}', 'bern_c5', 1e-2, 0., 0.1)

        pdfs['bern0'] = RooBernstein(f'bern0_{tag}', 'bern0', x, 
                                      RooArgList(params['bern_c0']))
        pdfs['bern1'] = RooBernstein(f'bern1_{tag}', 'bern1', x, 
                                      RooArgList(params['bern_c0'], params['bern_c1']))
        pdfs['bern2'] = RooBernstein(f'bern2_{tag}', 'bern2', x, 
                                      RooArgList(params['bern_c0'], params['bern_c1'], params['bern_c2']))
        pdfs['bern3'] = RooBernstein(f'bern3_{tag}', 'bern3', x, 
                                      RooArgList(params['bern_c0'], params['bern_c1'], 
                                                params['bern_c2'], params['bern_c3']))
        pdfs['bern4'] = RooBernstein(f'bern4_{tag}', 'bern4', x,
                                      RooArgList(params['bern_c0'], params['bern_c1'], 
                                                params['bern_c2'], params['bern_c3'], params['bern_c4']))
        pdfs['bern5'] = RooBernstein(f'bern5_{tag}', 'bern5', x,
                                      RooArgList(params['bern_c0'], params['bern_c1'], 
                                                params['bern_c2'], params['bern_c3'], 
                                                params['bern_c4'], params['bern_c5']))

        # ===== CHEBYCHEV PDFs =====
        params['chebychev_c0'] = RooRealVar(f'chebychev_c0_{tag}', 'chebychev_c0', 1.08, -1.1, 10.)
        params['chebychev_c1'] = RooRealVar(f'chebychev_c1_{tag}', 'chebychev_c1', 0.4, -1., 1.)
        params['chebychev_c2'] = RooRealVar(f'chebychev_c2_{tag}', 'chebychev_c2', 0.01, -0.1, 0.1)
        params['chebychev_c3'] = RooRealVar(f'chebychev_c3_{tag}', 'chebychev_c3', 0., -1., 1.)

        pdfs['chebychev1'] = RooChebychev(f'chebychev1_{tag}', 'chebychev1', x,
                                          RooArgList(params['chebychev_c0'], params['chebychev_c1']))
        pdfs['chebychev2'] = RooChebychev(f'chebychev2_{tag}', 'chebychev2', x,
                                          RooArgList(params['chebychev_c0'], params['chebychev_c1'], 
                                                    params['chebychev_c2']))
        pdfs['chebychev3'] = RooChebychev(f'chebychev3_{tag}', 'chebychev3', x,
                                          RooArgList(params['chebychev_c0'], params['chebychev_c1'], 
                                                    params['chebychev_c2'], params['chebychev_c3']))

        # ===== GAUSSIAN PDF =====
        params['gauss_mu'] = RooRealVar(f'gauss_mu{tag}', 'gauss_mu', 120, 100, 140)
        params['gauss_sigma'] = RooRealVar(f'gauss_sigma{tag}', 'gauss_sigma', 30, 20, 40)
        pdfs['gauss'] = RooGaussian(f'gauss{tag}', 'gauss', x, 
                                     params['gauss_mu'], params['gauss_sigma'])

        # ===== EXPONENTIAL PDFs =====
        params['exp_p1'] = RooRealVar(f'exp_p1{tag}', 'exp_p1', -0.0207, -0.022, -0.018)
        params['exp_p2'] = RooRealVar(f'exp_p2{tag}', 'exp_p2', -1e-2, -10, 10)
        params['exp_p3'] = RooRealVar(f'exp_p3{tag}', 'exp_p3', -1e-3, -10, 0)
        params['exp_c1'] = RooRealVar(f'exp_c1{tag}', 'exp_c1', 0., 1.)
        params['exp_c2'] = RooRealVar(f'exp_c2{tag}', 'exp_c2', 0., 1.)

        formulas['exp_frac1'] = RooFormulaVar(f"exp_frac1{tag}", "@0", 
                                              RooArgList(params['exp_c1']))
        formulas['exp_frac2'] = RooFormulaVar(f"exp_frac2{tag}", "(1-@0)*@1", 
                                              RooArgList(params['exp_c1'], params['exp_c2']))

        # Store intermediate exponential PDFs
        components['pdf_exp1'] = RooExponential(f'exp1_{tag}', 'exp1', x, params['exp_p1'])
        components['pdf_single_exp2'] = RooExponential(f'single_exp2{tag}', 'single_exp2', x, params['exp_p2'])
        components['pdf_single_exp3'] = RooExponential(f'single_exp3{tag}', 'single_exp3', x, params['exp_p3'])

        pdfs['exp1'] = components['pdf_exp1']
        pdfs['exp2'] = RooAddPdf(f'exp2{tag}', 'exp2',
                                 RooArgList(components['pdf_exp1'], components['pdf_single_exp2']),
                                 RooArgList(formulas['exp_frac1']))
        pdfs['exp3'] = RooAddPdf(f'exp3{tag}', 'exp3',
                                 RooArgList(components['pdf_exp1'], components['pdf_single_exp2'], 
                                           components['pdf_single_exp3']),
                                 RooArgList(formulas['exp_frac1'], formulas['exp_frac2']))

        # ===== POWER PDFs =====
        params['frac1'] = RooRealVar(f'frac1{tag}', 'frac1', 0.3, 0., 1.)
        params['frac2'] = RooRealVar(f'frac2{tag}', 'frac2', 0.5, 0., 1.)
        formulas['pow_frac1'] = RooFormulaVar(f"pow_frac1{tag}", "@0", 
                                              RooArgList(params['frac1']))
        formulas['pow_frac2'] = RooFormulaVar(f"pow_frac2{tag}", "(1-@0)*@1", 
                                              RooArgList(params['frac1'], params['frac2']))

        params['pow_p1'] = RooRealVar(f'p1{tag}', 'p1', -2.555, -10., 0.)
        params['pow_p2'] = RooRealVar(f'p2{tag}', 'p2', -8., -10., 0.)
        params['pow_p3'] = RooRealVar(f'p3{tag}', 'p3', -10., -10., 0.)

        # Store intermediate power PDFs
        components['pdf_pow1'] = RooGenericPdf(f'pow1_{tag}', 'pow1',
                                               'TMath::Power(@0,@1)',
                                               RooArgList(x, params['pow_p1']))
        components['pdf_pow2_comp'] = RooGenericPdf(f'pow2comp_{tag}', 'pow2comp',
                                                    'TMath::Power(@0,@1)',
                                                    RooArgList(x, params['pow_p2']))
        components['pdf_pow3_comp'] = RooGenericPdf(f'pow3comp_{tag}', 'pow3comp',
                                                    'TMath::Power(@0,@1)',
                                                    RooArgList(x, params['pow_p3']))

        pdfs['pow1'] = components['pdf_pow1']
        pdfs['pow2'] = RooAddPdf(f'pow2_{tag}', 'pow2',
                                 RooArgList(components['pdf_pow1'], components['pdf_pow2_comp']),
                                 RooArgList(formulas['pow_frac1']))
        pdfs['pow3'] = RooAddPdf(f'pow3_{tag}', 'pow3',
                                 RooArgList(components['pdf_pow1'], components['pdf_pow2_comp'], 
                                           components['pdf_pow3_comp']),
                                 RooArgList(formulas['pow_frac1'], formulas['pow_frac2']))

        # Return EVERYTHING to keep all references alive
        return {
            'pdfs': pdfs,
            'params': params,
            'components': components,
            'formulas': formulas
        }
