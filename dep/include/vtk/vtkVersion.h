/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkVersion.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkVersion
 * @brief   Versioning class for vtk
 *
 * Holds methods for defining/determining the current vtk version
 * (major, minor, build).
 *
 * @warning
 * This file will change frequently to update the VTKSourceVersion which
 * timestamps a particular source release.
 */

#ifndef vtkVersion_h
#define vtkVersion_h

#include "vtkCommonCoreModule.h" // For export macro
#include "vtkObject.h"
#include "vtkVersionMacros.h" // For version macros

#define VTK_SOURCE_VERSION "vtk version " VTK_VERSION

#define VTK_VERSION_CHECK(major, minor, build)                                                     \
  (10000000000ULL * major + 100000000ULL * minor + build)
#define VTK_VERSION_NUMBER                                                                         \
  VTK_VERSION_CHECK(VTK_MAJOR_VERSION, VTK_MINOR_VERSION, VTK_BUILD_VERSION)

class VTKCOMMONCORE_EXPORT vtkVersion : public vtkObject
{
public:
  static vtkVersion* New();
  vtkTypeMacro(vtkVersion, vtkObject);

  /**
   * Return the version of vtk this object is a part of.
   * A variety of methods are included. GetVTKSourceVersion returns a string
   * with an identifier which timestamps a particular source tree.
   */
  static const char* GetVTKVersion() { return VTK_VERSION; }
  static int GetVTKMajorVersion() { return VTK_MAJOR_VERSION; }
  static int GetVTKMinorVersion() { return VTK_MINOR_VERSION; }
  static int GetVTKBuildVersion() { return VTK_BUILD_VERSION; }
  static const char* GetVTKSourceVersion() { return VTK_SOURCE_VERSION; }

protected:
  vtkVersion() {} // insure constructor/destructor protected
  ~vtkVersion() override {}

private:
  vtkVersion(const vtkVersion&) = delete;
  void operator=(const vtkVersion&) = delete;
};

extern "C"
{
  VTKCOMMONCORE_EXPORT const char* GetVTKVersion();
}

#endif

// VTK-HeaderTest-Exclude: vtkVersion.h
