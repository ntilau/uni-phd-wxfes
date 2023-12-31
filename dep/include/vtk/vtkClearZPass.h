/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkClearZPass.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkClearZPass
 * @brief   Clear the depth buffer with a given value.
 *
 * Clear the depth buffer with a given value.
 *
 * @sa
 * vtkRenderPass
 */

#ifndef vtkClearZPass_h
#define vtkClearZPass_h

#include "vtkRenderPass.h"
#include "vtkRenderingOpenGL2Module.h" // For export macro

class vtkOpenGLRenderWindow;

class VTKRENDERINGOPENGL2_EXPORT vtkClearZPass : public vtkRenderPass
{
public:
  static vtkClearZPass* New();
  vtkTypeMacro(vtkClearZPass, vtkRenderPass);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  /**
   * Perform rendering according to a render state \p s.
   * \pre s_exists: s!=0
   */
  void Render(const vtkRenderState* s) override;

  //@{
  /**
   * Set/Get the depth value. Initial value is 1.0 (farest).
   */
  vtkSetClampMacro(Depth, double, 0.0, 1.0);
  vtkGetMacro(Depth, double);
  //@}

protected:
  /**
   * Default constructor.
   */
  vtkClearZPass();

  /**
   * Destructor.
   */
  ~vtkClearZPass() override;

  double Depth;

private:
  vtkClearZPass(const vtkClearZPass&) = delete;
  void operator=(const vtkClearZPass&) = delete;
};

#endif
