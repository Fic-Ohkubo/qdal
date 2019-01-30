/* external C library should implement the following functions */
extern QErrorCode Matrix_Create(LANG_C_MATRIX*);
extern QErrorCode Matrix_Destroy(LANG_C_MATRIX*);
#if defined(ADAPTER_BLOCK_REAL)
extern QErrorCode Matrix_BlockCreate(LANG_C_MATRIX*,const QInt);
#endif
extern QErrorCode Matrix_SetSymType(LANG_C_MATRIX*,const QcSymType);
#if defined(ADAPTER_BLOCK_REAL)
extern QErrorCode Matrix_SetNonZeroBlocks(LANG_C_MATRIX*,
                                          const QInt,
                                          const QInt[],
                                          const QInt[]);
#elif defined(ADAPTER_CMPLX_MAT)
extern QErrorCode Matrix_SetDataType(LANG_C_MATRIX*,const QcDataType);
#endif
extern QErrorCode Matrix_SetDimMat(LANG_C_MATRIX*,const QInt);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode Matrix_SetStorageMode(LANG_C_MATRIX*,const QcStorageMode);
#endif
extern QErrorCode Matrix_Assemble(LANG_C_MATRIX*);
#if defined(ADAPTER_BLOCK_REAL)
extern QErrorCode Matrix_GetDimBlock(LANG_C_MATRIX*,QInt*);
#endif
extern QErrorCode Matrix_GetSymType(LANG_C_MATRIX*,QcSymType*);
#if defined(ADAPTER_BLOCK_REAL)
extern QErrorCode Matrix_GetNonZeroBlocks(LANG_C_MATRIX*,
                                          const QInt,
                                          const QInt[],
                                          const QInt[],
                                          QBool*);
#elif defined(ADAPTER_CMPLX_MAT)
extern QErrorCode Matrix_GetDataType(LANG_C_MATRIX*,QcDataType*);
#endif
extern QErrorCode Matrix_GetDimMat(LANG_C_MATRIX*,QInt*);
#if defined(QCMATRIX_STORAGE_MODE)
extern QErrorCode Matrix_GetStorageMode(LANG_C_MATRIX*,QcStorageMode*);
#endif
extern QErrorCode Matrix_IsAssembled(LANG_C_MATRIX*,QBool*);
#if defined(ADAPTER_BLOCK_REAL)
extern QErrorCode Matrix_SetValues(LANG_C_MATRIX*,
                                   const QInt,
                                   const QInt,
                                   const QInt[],
                                   const QInt[],
                                   const QReal*);
extern QErrorCode Matrix_GetValues(LANG_C_MATRIX*,
                                   const QInt,
                                   const QInt,
                                   const QInt[],
                                   const QInt[],
                                   QReal*);
#elif defined(ADAPTER_CMPLX_MAT)
extern QErrorCode Matrix_SetValues(LANG_C_MATRIX*,
                                   const QInt[],
                                   const QInt[],
                                   const QReal*,
                                   const QReal*);
extern QErrorCode Matrix_GetValues(LANG_C_MATRIX*,
                                   const QInt[],
                                   const QInt[],
                                   QReal*,
                                   QReal*);
#elif defined(ADAPTER_REAL_MAT)
extern QErrorCode Matrix_SetValues(LANG_C_MATRIX*,
                                   const QInt[],
                                   const QInt[],
                                   const QReal*);
extern QErrorCode Matrix_GetValues(LANG_C_MATRIX*,
                                   const QInt[],
                                   const QInt[],
                                   QReal*);
#endif
extern QErrorCode Matrix_Duplicate(LANG_C_MATRIX*,const QcDuplicateOption,LANG_C_MATRIX*);
extern QErrorCode Matrix_ZeroEntries(LANG_C_MATRIX*);
#if defined(ADAPTER_BLOCK_REAL)
extern QErrorCode Matrix_GetTrace(LANG_C_MATRIX*,const QInt,QReal*);
#elif defined(ADAPTER_CMPLX_MAT)
extern QErrorCode Matrix_GetTrace(LANG_C_MATRIX*,QReal*);
#elif defined(ADAPTER_REAL_MAT)
extern QErrorCode Matrix_GetTrace(LANG_C_MATRIX*,QReal*);
#endif
#if defined(QCMATRIX_ENABLE_VIEW)
extern QErrorCode Matrix_Write(LANG_C_MATRIX*,const QChar*,const QcViewOption);
extern QErrorCode Matrix_Read(LANG_C_MATRIX*,const QChar*,const QcViewOption);
#endif
#if defined(ADAPTER_CMPLX_MAT)
extern QErrorCode Matrix_Scale(const QReal[],LANG_C_MATRIX*);
extern QErrorCode Matrix_AXPY(const QReal[],LANG_C_MATRIX*,LANG_C_MATRIX*);
#elif defined(ADAPTER_BLOCK_REAL) || defined(ADAPTER_REAL_MAT)
extern QErrorCode Matrix_Scale(const QReal,LANG_C_MATRIX*);
extern QErrorCode Matrix_AXPY(const QReal,LANG_C_MATRIX*,LANG_C_MATRIX*);
#endif
extern QErrorCode Matrix_Transpose(const QcMatOperation,LANG_C_MATRIX*,LANG_C_MATRIX*);
#if defined(ADAPTER_CMPLX_MAT)
extern QErrorCode Matrix_GEMM(const QcMatOperation,
                              const QcMatOperation,
                              const QReal[],
                              LANG_C_MATRIX*,
                              LANG_C_MATRIX*,
                              const QReal[],
                              LANG_C_MATRIX*);
#elif defined(ADAPTER_BLOCK_REAL) || defined(ADAPTER_REAL_MAT)
extern QErrorCode Matrix_GEMM(const QcMatOperation,
                              const QcMatOperation,
                              const QReal,
                              LANG_C_MATRIX*,
                              LANG_C_MATRIX*,
                              const QReal,
                              LANG_C_MATRIX*);
#endif
