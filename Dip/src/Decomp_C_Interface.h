#ifndef DECOMP_C_INTERFACE_INCLUDED
#define DECOMP_C_INTERFACE_INCLUDED

#ifdef __cplusplus
extern "C"
{
#endif
    typedef struct UtilParameters UtilParameters;
    typedef struct DecompApp DecompApp;
    typedef struct DecompConstraintSet DecompConstraintSet;
    typedef struct DecompAlgo DecompAlgo;
    typedef struct AlpsDecompModel AlpsDecompModel;
    typedef struct DecompVar DecompVar;

    //===========================================================================//
    // UtilParameters
    //===========================================================================//
    UtilParameters *Dip_UtilParameters_new();
    void Dip_UtilParameters_delete(UtilParameters *utilParams);
    void Dip_UtilParameters_ScanCmdLineArgs(UtilParameters *utilParams,
                                            int argc,
                                            char *argv[]);
    int Dip_UtilParameters_getSetting(UtilParameters *utilParams,
                                      const char *name,
                                      const int defaultValue,
                                      const char *section);

    //===========================================================================//
    // DecompApp
    //===========================================================================//
    DecompApp *Dip_DecompApp_new(UtilParameters *utilParams);
    void Dip_DecompApp_delete(DecompApp *app);
    void Dip_DecompApp_setModelObjective(DecompApp *app, const double *objective, const int length);
    void Dip_DecompApp_setModelCore(DecompApp *app, DecompConstraintSet *model,
                                    const char *modelName);
    void Dip_DecompApp_setModelRelax(DecompApp *app,
                                     DecompConstraintSet *model,
                                     const char *modelName,
                                     const int blockId);
    double Dip_DecompApp_getInfinity(DecompApp *app);

    //===========================================================================//
    // DecompConstraintSet
    //===========================================================================//
    DecompConstraintSet *Dip_DecompConstraintSet_new();
    void Dip_DecompConstraintSet_delete(DecompConstraintSet *constrSet);
    void Dip_DecompConstraintSet_init(DecompConstraintSet *constrSet, const int nCols, const int nRows);
    void Dip_DecompConstraintSet_appendRow(DecompConstraintSet *constrSet, int size, const int *inds, const double *elems,
                                           double loBound,
                                           double upBound,
                                           const char *rowName);
    void Dip_DecompConstraintSet_pushCol(DecompConstraintSet *constrSet, const double loBound,
                                         const double upBound,
                                         const char *colName,
                                         const int isInteger,
                                         const int origIndex);
    void Dip_DecompConstraintSet_setSparse(DecompConstraintSet *constrSet, const int numColsOrig);

    //===========================================================================//
    // DecompAlgo
    //===========================================================================//
    DecompAlgo *Dip_DecompAlgoC_new(DecompApp *app, UtilParameters *utilParam);
    DecompAlgo *Dip_DecompAlgoPC_new(DecompApp *app, UtilParameters *utilParam);
    void Dip_DecompAlgo_delete(DecompAlgo *algo);

    //===========================================================================//
    // AlpsDecompModel
    //===========================================================================//
    AlpsDecompModel *Dip_AlpsDecompModel_new(UtilParameters *utilParam, DecompAlgo *algo);
    void Dip_AlpsDecompModel_delete(AlpsDecompModel *alps);
    int Dip_AlpsDecompModel_solve(AlpsDecompModel *alps);

    //===========================================================================//
    // DecompVar
    //===========================================================================//
    DecompVar *Dip_DecompVar_new(const int len,
                                 const int *ind,
                                 const double *els,
                                 const double redCost,
                                 const double origCost,
                                 const int varType);
    void Dip_DecompVar_delete(DecompVar *var);

#ifdef __cplusplus
}
#endif
#endif
