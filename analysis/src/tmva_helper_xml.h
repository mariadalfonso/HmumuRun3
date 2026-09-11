#ifndef tmva_xml_h
#define tmva_xml_h

#include "ROOT/RVec.hxx"

using Vec_f = ROOT::VecOps::RVec<float>;

class tmva_xml {

    public:
        tmva_xml(const std::string &filename) {

            auto c = TMVA::Experimental::Internal::ParseXMLConfig(filename);
			fVariables = c.variables;
			fExpressions = c.variable_expressions;
			fAnalysisType = c.analysisType;
			fNumClasses = c.numClasses;

            fReader = std::make_unique<TMVA::Reader>("Silent");
			const auto numVars = fVariables.size();
            fValues = std::vector<float>(numVars);

            for(std::size_t i = 0; i < numVars; i++) {
                fReader->AddVariable(TString(fExpressions[i]), &fValues[i]);
            }
            fReader->BookMVA(name, filename.c_str());

        }
        ~tmva_xml() {};
  
        std::vector<float> Compute(const Vec_f &x) {

            if (x.size() != fVariables.size())
                throw std::runtime_error("Size of input vector is not equal to number of variables.");

            // Copy over inputs to memory used by TMVA reader
            for (std::size_t i = 0; i < x.size(); i++) {
                fValues[i] = x[i];
            }

            // Evaluate TMVA model
            // Classification
            if (fAnalysisType == TMVA::Experimental::Internal::AnalysisType::Classification) {
                return std::vector<float>({static_cast<float>(fReader->EvaluateMVA(name))});
            }
            // Regression
            else if (fAnalysisType == TMVA::Experimental::Internal::AnalysisType::Regression) {
                return fReader->EvaluateRegression(name);
            }
            // Multiclass
            else if (fAnalysisType == TMVA::Experimental::Internal::AnalysisType::Multiclass) {
                return fReader->EvaluateMulticlass(name);
            }
            // Throw error
            else {
                throw std::runtime_error("RReader has undefined analysis type.");
                return std::vector<float>();
            }
        }

    private:
        std::unique_ptr<TMVA::Reader> fReader;
		std::vector<float> fValues;
		std::vector<std::string> fVariables;
		std::vector<std::string> fExpressions;
		unsigned int fNumClasses;
		const char *name = "RReader";
		TMVA::Experimental::Internal::AnalysisType fAnalysisType;
};

#endif

#ifndef tmva_helper_xml_h
#define tmva_helper_xml_h

#include "ROOT/RVec.hxx"

using Vec_f = ROOT::VecOps::RVec<float>;

class tmva_helper_xml {
    public:
        tmva_helper_xml(const std::string &filename, const unsigned int nslots = 1) {

            const unsigned int nslots_actual = std::max(nslots, 1U);

			for (unsigned int islot = 0; islot < nslots_actual; ++islot) {
			  interpreters_.emplace_back(std::make_shared<tmva_xml>(filename));
            }
        }

    ~tmva_helper_xml() {};
    std::vector<float> operator()(unsigned int slot, const Vec_f & vars) {
 			return interpreters_[slot]->Compute(vars);
        }

    private:
  std::vector<std::shared_ptr<tmva_xml>> interpreters_;
};


#endif


#ifndef tmva_helper_xml_vec_h
#define tmva_helper_xml_vec_h

#include "ROOT/RVec.hxx"

using Vec_f = ROOT::VecOps::RVec<float>;
using Vec_Vec_f = ROOT::VecOps::RVec<Vec_f>;

template <typename... Args>
Vec_Vec_f makeTMVAVarsVec(const Args&... args)
{
    auto inputs = std::forward_as_tuple(args...);

    Vec_Vec_f out;

    if constexpr (sizeof...(Args) == 0)
        return out;

    const size_t nlep = std::get<0>(inputs).size();
    out.reserve(nlep);

    for (size_t i = 0; i < nlep; ++i) {

        Vec_f vars;
        vars.reserve(sizeof...(Args));

        std::apply([&](const auto&... v) {
            (vars.push_back(v[i]), ...);
        }, inputs);

        out.push_back(std::move(vars));
    }

    return out;
}

class tmva_helper_xml_vec {
public:

  tmva_helper_xml_vec(const std::string &filename,
		      const unsigned int nslots = 1)
  {
    const unsigned int nslots_actual = std::max(nslots, 1U);

    for (unsigned int islot = 0; islot < nslots_actual; ++islot) {
      interpreters_.emplace_back(std::make_shared<tmva_xml>(filename));
    }
  }


  Vec_f operator()(
		   unsigned int slot,
		   const Vec_Vec_f& vars
		   )
  {
    Vec_f scores(vars.size());

    for (size_t i = 0; i < vars.size(); ++i) {
      auto result = interpreters_[slot]->Compute(vars[i]);
      scores[i] = result[0];
    }

    return scores;
  }

private:
  std::vector<std::shared_ptr<tmva_xml>> interpreters_;
};

#endif
