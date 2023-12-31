/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkXMLPartitionedDataSetWriter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkXMLPartitionedDataSetWriter
 * @brief   writer for vtkPartitionedDataSet.
 *
 * vtkXMLPartitionedDataSetWriter is a vtkXMLCompositeDataWriter subclass to handle
 * vtkPartitionedDataSet.
 */

#ifndef vtkXMLPartitionedDataSetWriter_h
#define vtkXMLPartitionedDataSetWriter_h

#include "vtkIOXMLModule.h" // For export macro
#include "vtkXMLCompositeDataWriter.h"

class VTKIOXML_EXPORT vtkXMLPartitionedDataSetWriter : public vtkXMLCompositeDataWriter
{
public:
  static vtkXMLPartitionedDataSetWriter* New();
  vtkTypeMacro(vtkXMLPartitionedDataSetWriter, vtkXMLCompositeDataWriter);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /**
   * Get the default file extension for files written by this writer.
   */
  const char* GetDefaultFileExtension() override { return "vtpd"; }

protected:
  vtkXMLPartitionedDataSetWriter();
  ~vtkXMLPartitionedDataSetWriter() override;

  int FillInputPortInformation(int port, vtkInformation* info) override;

  // Internal method called recursively to create the xml tree for the children
  // of compositeData.
  int WriteComposite(
    vtkCompositeDataSet* compositeData, vtkXMLDataElement* parent, int& writerIdx) override;

private:
  vtkXMLPartitionedDataSetWriter(const vtkXMLPartitionedDataSetWriter&) = delete;
  void operator=(const vtkXMLPartitionedDataSetWriter&) = delete;
};

#endif
