#ifdef _MSC_VER
  #if defined(BUILD_DLL) || defined(USE_DLL)
    #ifdef BUILD_DLL
      #define PEL_EXPORT __declspec(dllexport)
    #else
      #define PEL_EXPORT __declspec(dllimport)
    #endif
  #else
    #define PEL_EXPORT
  #endif
#else
  #define PEL_EXPORT
#endif
