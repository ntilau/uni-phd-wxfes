/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkExtractBlock.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkExtractBlock
 * @brief   extracts blocks from a multiblock dataset.
 *
 * vtkExtractBlock is a filter that extracts blocks from a multiblock
 * dataset.  Each node in the multi-block tree is identified by an \c
 * index. The index can be obtained by performing a preorder traversal of the
 * tree (including empty nodes). eg. A(B (D, E), C(F, G)).  Inorder traversal
 * yields: A, B, D, E, C, F, G Index of A is 0, while index of C is 4.
 *
 * Note that if you specify node 0, then the input is simply shallow copied
 * to the output. This is true even if other nodes are specified along with
 * node 0.
 */

#ifndef vtkExtractBlock_h
#define vtkExtractBlock_h

#include "vtkFiltersExtractionModule.h" // For export macro
#include "vtkMultiBlockDataSetAlgorithm.h"

class vtkDataObjectTreeIterator;
class vtkMultiPieceDataSet;

class VTKFILTERSEXTRACTION_EXPORT vtkExtractBlock : public vtkMultiBlockDataSetAlgorithm
{
public:
  //@{
  /**
   * Standard methods for instantiation, type information, and printing.
   */
  static vtkExtractBlock* New();
  vtkTypeMacro(vtkExtractBlock, vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) override;
  //@{

  //@{
  /**
   * Select the block indices to extract.  Each node in the multi-block tree
   * is identified by an \c index. The index can be obtained by performing a
   * preorder traversal of the tree (including empty nodes). eg. A(B (D, E),
   * C(F, G)).  Inorder traversal yields: A, B, D, E, C, F, G Index of A is
   * 0, while index of C is 4. (Note: specifying node 0 means the input is
   * copied to the output.)
   */
  void AddIndex(unsigned int index);
  void RemoveIndex(unsigned int index);
  void RemoveAllIndices();
  //@}

  //@{
  /**
   * When set, the output multiblock dataset will be pruned to remove empty
   * nodes. On by default.
   */
  vtkSetMacro(PruneOutput, vtkTypeBool);
  vtkGetMacro(PruneOutput, vtkTypeBool);
  vtkBooleanMacro(PruneOutput, vtkTypeBool);
  //@}

  //@{
  /**
   * This is used only when PruneOutput is ON. By default, when pruning the
   * output i.e. remove empty blocks, if node has only 1 non-null child block,
   * then that node is removed. To preserve these parent nodes, set this flag to
   * true. Off by default.
   */
  vtkSetMacro(MaintainStructure, vtkTypeBool);
  vtkGetMacro(MaintainStructure, vtkTypeBool);
  vtkBooleanMacro(MaintainStructure, vtkTypeBool);
  //@}

protected:
  vtkExtractBlock();
  ~vtkExtractBlock() override;

  /**
   * Internal key, used to avoid pruning of a branch.
   */
  static vtkInformationIntegerKey* DONT_PRUNE();

  /// Implementation of the algorithm.
  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*) override;

  /// Extract subtree
  void CopySubTree(
    vtkDataObjectTreeIterator* loc, vtkMultiBlockDataSet* output, vtkMultiBlockDataSet* input);
  bool Prune(vtkMultiBlockDataSet* mblock);
  bool Prune(vtkMultiPieceDataSet* mblock);
  bool Prune(vtkDataObject* mblock);

  vtkTypeBool PruneOutput;
  vtkTypeBool MaintainStructure;

private:
  vtkExtractBlock(const vtkExtractBlock&) = delete;
  void operator=(const vtkExtractBlock&) = delete;

  class vtkSet;
  vtkSet* Indices;
  vtkSet* ActiveIndices;
};

#endif
