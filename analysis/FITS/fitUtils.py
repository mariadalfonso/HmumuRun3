# fitUtils.py
import ROOT
from ROOT import RooFit, RooRealVar, RooDataHist, RooArgList

class FitUtilities:
    """Centralized fitting utilities"""
    
    @staticmethod
    def create_fit_cmd_list(minimizer="Minuit2", strategy=2):
        """Create standard RooFit fit command list"""
        cmdList = ROOT.RooLinkedList()
        cmdList.Add(RooFit.Minimizer(minimizer))
        cmdList.Add(RooFit.Strategy(strategy))
        cmdList.Add(RooFit.Range("full"))
        cmdList.Add(RooFit.Save(True))
        return cmdList

    @staticmethod
    def fit_model(model, data, cmdList=None, data_blinded=None):
        """
        Fit a model to data
        
        Args:
            model: RooAbsPdf to fit
            data: RooDataHist or RooDataSet
            cmdList: RooLinkedList of fit options
            data_blinded: optional blinded data
            
        Returns:
            RooFitResult
        """
        if cmdList is None:
            cmdList = FitUtilities.create_fit_cmd_list()
        
        fit_data = data_blinded if data_blinded is not None else data
        fitresult = model.fitTo(fit_data, cmdList)
        return fitresult

    @staticmethod
    def calculate_chi2_from_fit(fitresult, n_data_points):
        """
        Calculate chi2 and chi2/ndof from fit result
        
        Args:
            fitresult: RooFitResult
            n_data_points: number of data points
            
        Returns:
            dict with 'chi2', 'ndof', 'chi2_ndof'
        """
        chi2 = fitresult.minNll() * 2  # Convert NLL to chi2
        n_params = fitresult.floatParsFinal().getSize()
        ndof = n_data_points - n_params
        chi2_ndof = chi2 / ndof if ndof > 0 else -1.0
        
        return {
            'chi2': chi2,
            'ndof': ndof,
            'chi2_ndof': chi2_ndof,
            'n_params': n_params
        }

    @staticmethod
    def print_fit_results(model_name, fitresult, n_data_points, file_obj=None):
        """
        Print fit results in standard format
        
        Args:
            model_name: name of model
            fitresult: RooFitResult
            n_data_points: number of data points
            file_obj: optional file to write to
        """
        chi2_info = FitUtilities.calculate_chi2_from_fit(fitresult, n_data_points)
        
        output_str = (f"{model_name:30s}  chi2/ndof={chi2_info['chi2_ndof']:7.3f}  "
                     f"ndof={chi2_info['ndof']:3d}  n_params={chi2_info['n_params']:2d}\n")
        
        print(output_str.rstrip('\n'))
        if file_obj:
            file_obj.write(output_str)
        
        return chi2_info

    @staticmethod
    def select_best_model(fitresults_dict):
        """
        Select best model based on lowest chi2/ndof
        
        Args:
            fitresults_dict: dict with model names as keys and chi2_info as values
            
        Returns:
            tuple (best_model_name, chi2_info)
        """
        best_model = min(fitresults_dict.items(), 
                        key=lambda x: x[1]['chi2_ndof'] if x[1]['chi2_ndof'] > 0 else float('inf'))
        return best_model
