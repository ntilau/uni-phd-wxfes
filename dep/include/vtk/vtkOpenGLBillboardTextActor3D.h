/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOpenGLBillboardTextActor3D.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class vtkOpenGLBillboardTextActor3D
 * @brief Handles GL2PS capture of billboard text.
 */

#ifndef vtkOpenGLBillboardTextActor3D_h
#define vtkOpenGLBillboardTextActor3D_h

#include "vtkBillboardTextActor3D.h"
#include "vtkRenderingOpenGL2Module.h" // For export macro

class vtkOpenGLGL2PSHelper;

class VTKRENDERINGOPENGL2_EXPORT vtkOpenGLBillboardTextActor3D : public vtkBillboardTextActor3D
{
public:
  static vtkOpenGLBillboardTextActor3D* New();
  vtkTypeMacro(vtkOpenGLBillboardTextActor3D, vtkBillboardTextActor3D);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  int RenderTranslucentPolygonalGeometry(vtkViewport* vp) override;

protected:
  vtkOpenGLBillboardTextActor3D();
  ~vtkOpenGLBillboardTextActor3D() override;

  int RenderGL2PS(vtkViewport* viewport, vtkOpenGLGL2PSHelper* gl2ps);

private:
  vtkOpenGLBillboardTextActor3D(const vtkOpenGLBillboardTextActor3D&) = delete;
  void operator=(const vtkOpenGLBillboardTextActor3D&) = delete;
};

#endif // vtkOpenGLBillboardTextActor3D_h
