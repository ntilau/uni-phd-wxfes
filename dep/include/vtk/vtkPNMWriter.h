/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPNMWriter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkPNMWriter
 * @brief   Writes PNM (portable any map)  files.
 *
 * vtkPNMWriter writes PNM file. The data type
 * of the file is unsigned char regardless of the input type.
 */

#ifndef vtkPNMWriter_h
#define vtkPNMWriter_h

#include "vtkIOImageModule.h" // For export macro
#include "vtkImageWriter.h"

class VTKIOIMAGE_EXPORT vtkPNMWriter : public vtkImageWriter
{
public:
  static vtkPNMWriter* New();
  vtkTypeMacro(vtkPNMWriter, vtkImageWriter);
  void PrintSelf(ostream& os, vtkIndent indent) override;

protected:
  vtkPNMWriter() {}
  ~vtkPNMWriter() override {}

  void WriteFile(ostream* file, vtkImageData* data, int extent[6], int wExt[6]) override;
  void WriteFileHeader(ostream*, vtkImageData*, int wExt[6]) override;

private:
  vtkPNMWriter(const vtkPNMWriter&) = delete;
  void operator=(const vtkPNMWriter&) = delete;
};

#endif
