/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkLogLookupTable.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkLogLookupTable
 * @brief   map scalars into colors using log (base 10) scale
 *
 * This class is an empty shell.  Use vtkLookupTable with SetScaleToLog10()
 * instead.
 *
 * @sa
 * vtkLookupTable
 */

#ifndef vtkLogLookupTable_h
#define vtkLogLookupTable_h

#include "vtkLookupTable.h"
#include "vtkRenderingCoreModule.h" // For export macro

class VTKRENDERINGCORE_EXPORT vtkLogLookupTable : public vtkLookupTable
{
public:
  static vtkLogLookupTable* New();

  vtkTypeMacro(vtkLogLookupTable, vtkLookupTable);
  void PrintSelf(ostream& os, vtkIndent indent) override;

protected:
  vtkLogLookupTable(int sze = 256, int ext = 256);
  ~vtkLogLookupTable() override {}

private:
  vtkLogLookupTable(const vtkLogLookupTable&) = delete;
  void operator=(const vtkLogLookupTable&) = delete;
};

#endif
