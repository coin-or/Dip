#include "Decomp_C_Interface.h"
#include "DecompApp.h"
#include "AlpsDecompModel.h"
#include "DecompAlgoC.h"
#include "DecompAlgoPC.h"

extern "C"
{
    //===========================================================================//
    // UtilParameters
    //===========================================================================//
    UtilParameters *Dip_UtilParameters_new()
    {
        return new UtilParameters();
    }

    void Dip_UtilParameters_delete(UtilParameters *utilParams)
    {
        if (utilParams)
            delete utilParams;
    }

    void Dip_UtilParameters_ScanCmdLineArgs(UtilParameters *utilParams,
                                            int argc,
                                            char *argv[])
    {
        utilParams->ScanCmdLineArgs(argc, argv);
    }

    int Dip_UtilParameters_getSetting(UtilParameters *utilParams,
                                      const char *name,
                                      const int defaultValue,
                                      const char *section)
    {
        return utilParams->GetSetting(name, defaultValue, section);
    }

    //===========================================================================//
    // DecompApp
    //===========================================================================//
    DecompApp *Dip_DecompApp_new(UtilParameters *utilParams)
    {
        return new DecompApp(*utilParams);
    }

    void Dip_DecompApp_delete(DecompApp *app)
    {
        if (app)
            delete app;
    }

    void Dip_DecompApp_setModelObjective(DecompApp *app, const double *objective, const int length)
    {
        app->setModelObjective(objective, length);
    }

    void Dip_DecompApp_setModelCore(DecompApp *app, DecompConstraintSet *model,
                                    const char *modelName)
    {
        std::string s(modelName);
        app->setModelCore(model, s);
    }

    void Dip_DecompApp_setModelRelax(DecompApp *app,
                                     DecompConstraintSet *model,
                                     const char *modelName,
                                     const int blockId)
    {
        std::string s(modelName);
        app->setModelRelax(model, s, blockId);
    }

    double Dip_DecompApp_getInfinity(DecompApp *app)
    {
        return app->m_infinity;
    }

    //===========================================================================//
    // DecompConstraintSet
    //===========================================================================//
    DecompConstraintSet *Dip_DecompConstraintSet_new()
    {
        return new DecompConstraintSet();
    }
    void Dip_DecompConstraintSet_delete(DecompConstraintSet *constrSet)
    {
        if (constrSet)
            delete constrSet;
    }
    void Dip_DecompConstraintSet_init(DecompConstraintSet *constrSet, const int nCols, const int nRows)
    {
        constrSet->init(nCols, nRows);
    }
    void Dip_DecompConstraintSet_appendRow(DecompConstraintSet *constrSet, int size, const int *inds, const double *elems,
                                           double loBound,
                                           double upBound,
                                           const char *rowName)
    {
        CoinPackedVector row(size, inds, elems);
        std::string s(rowName);
        constrSet->appendRow(row, loBound, upBound, s);
    }
    void Dip_DecompConstraintSet_pushCol(DecompConstraintSet *constrSet, const double loBound,
                                         const double upBound,
                                         const char *colName,
                                         const int isInteger,
                                         const int origIndex)
    {
        std::string s(colName);
        constrSet->pushCol(loBound, upBound, s, isInteger, origIndex);
    }
    void Dip_DecompConstraintSet_setSparse(DecompConstraintSet *constrSet, const int numColsOrig)
    {
        constrSet->setSparse(numColsOrig);
    }

    //===========================================================================//
    // DecompAlgo
    //===========================================================================//
    DecompAlgo *Dip_DecompAlgoC_new(DecompApp *app, UtilParameters *utilParam)
    {
        return new DecompAlgoC(app, *utilParam);
    }
    DecompAlgo *Dip_DecompAlgoPC_new(DecompApp *app, UtilParameters *utilParam)
    {
        return new DecompAlgoPC(app, *utilParam);
    }
    void Dip_DecompAlgo_delete(DecompAlgo *algo)
    {
        if (algo)
            delete algo;
    }

    //===========================================================================//
    // AlpsDecompModel
    //===========================================================================//
    AlpsDecompModel *Dip_AlpsDecompModel_new(UtilParameters *utilParam, DecompAlgo *algo)
    {
        return new AlpsDecompModel(*utilParam, algo);
    }
    void Dip_AlpsDecompModel_delete(AlpsDecompModel *alps)
    {
        if (alps)
            delete alps;
    }
    int Dip_AlpsDecompModel_solve(AlpsDecompModel *alps)
    {
        FILE *fp = fopen("master.lp", "w");
        alps->getDecompAlgo()->getMasterOSI()->writeLp(fp);
        AlpsExitStatus status = alps->solve();
        return status;
    }

    //===========================================================================//
    // DecompVar
    //===========================================================================//
    DecompVar *Dip_DecompVar_new(const int len,
                                 const int *ind,
                                 const double *els,
                                 const double redCost,
                                 const double origCost,
                                 const int varType)
    {
        return new DecompVar(len, ind, els, redCost, origCost, (DecompVarType)varType);
    }

    void Dip_DecompVar_delete(DecompVar *var)
    {
        delete var;
    }
}