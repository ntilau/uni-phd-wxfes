/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSMPWarpVector.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkSMPWarpVector
 * @brief   multithreaded vtkWarpVector
 *
 * Just like parent, but uses the SMP framework to do the work on many threads.
 */

#ifndef vtkSMPWarpVector_h
#define vtkSMPWarpVector_h

#include "vtkFiltersSMPModule.h" // For export macro
#include "vtkWarpVector.h"

class vtkInformation;
class vtkInformationVector;

#if !defined(VTK_LEGACY_REMOVE)
class VTKFILTERSSMP_EXPORT vtkSMPWarpVector : public vtkWarpVector
{
public:
  vtkTypeMacro(vtkSMPWarpVector, vtkWarpVector);
  static vtkSMPWarpVector* New();
  void PrintSelf(ostream& os, vtkIndent indent) override;

protected:
  vtkSMPWarpVector();
  ~vtkSMPWarpVector() override;

  /**
   * Overridden to use threads.
   */
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

private:
  vtkSMPWarpVector(const vtkSMPWarpVector&) = delete;
  void operator=(const vtkSMPWarpVector&) = delete;
};

#endif // VTK_LEGACY_REMOVE
#endif // vtkSMPWarpVector_h
