/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkDataSetWriter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkDataSetWriter
 * @brief   write any type of vtk dataset to file
 *
 * vtkDataSetWriter is an abstract class for mapper objects that write their
 * data to disk (or into a communications port). The input to this object is
 * a dataset of any type.
 */

#ifndef vtkDataSetWriter_h
#define vtkDataSetWriter_h

#include "vtkDataWriter.h"
#include "vtkIOLegacyModule.h" // For export macro

class VTKIOLEGACY_EXPORT vtkDataSetWriter : public vtkDataWriter
{
public:
  static vtkDataSetWriter* New();
  vtkTypeMacro(vtkDataSetWriter, vtkDataWriter);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  //@{
  /**
   * Get the input to this writer.
   */
  vtkDataSet* GetInput();
  vtkDataSet* GetInput(int port);
  //@}

protected:
  vtkDataSetWriter() {}
  ~vtkDataSetWriter() override {}

  void WriteData() override;

  int FillInputPortInformation(int port, vtkInformation* info) override;

private:
  vtkDataSetWriter(const vtkDataSetWriter&) = delete;
  void operator=(const vtkDataSetWriter&) = delete;
};

#endif
