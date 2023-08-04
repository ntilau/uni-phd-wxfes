
#ifndef VTKIOLSDYNA_EXPORT_H
#define VTKIOLSDYNA_EXPORT_H

#ifdef VTKIOLSDYNA_STATIC_DEFINE
#  define VTKIOLSDYNA_EXPORT
#  define VTKIOLSDYNA_NO_EXPORT
#else
#  ifndef VTKIOLSDYNA_EXPORT
#    ifdef IOLSDyna_EXPORTS
        /* We are building this library */
#      define VTKIOLSDYNA_EXPORT 
#    else
        /* We are using this library */
#      define VTKIOLSDYNA_EXPORT 
#    endif
#  endif

#  ifndef VTKIOLSDYNA_NO_EXPORT
#    define VTKIOLSDYNA_NO_EXPORT 
#  endif
#endif

#ifndef VTKIOLSDYNA_DEPRECATED
#  define VTKIOLSDYNA_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef VTKIOLSDYNA_DEPRECATED_EXPORT
#  define VTKIOLSDYNA_DEPRECATED_EXPORT VTKIOLSDYNA_EXPORT VTKIOLSDYNA_DEPRECATED
#endif

#ifndef VTKIOLSDYNA_DEPRECATED_NO_EXPORT
#  define VTKIOLSDYNA_DEPRECATED_NO_EXPORT VTKIOLSDYNA_NO_EXPORT VTKIOLSDYNA_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef VTKIOLSDYNA_NO_DEPRECATED
#    define VTKIOLSDYNA_NO_DEPRECATED
#  endif
#endif

#endif /* VTKIOLSDYNA_EXPORT_H */
