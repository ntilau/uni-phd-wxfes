/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkOpenGLTextMapper.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

/**
 * @class   vtkOpenGLTextMapper
 * @brief   vtkTextMapper override for OpenGL2.
 */

#ifndef vtkOpenGLTextMapper_h
#define vtkOpenGLTextMapper_h

#include "vtkRenderingOpenGL2Module.h" // For export macro
#include "vtkTextMapper.h"

class vtkOpenGLGL2PSHelper;

class VTKRENDERINGOPENGL2_EXPORT vtkOpenGLTextMapper : public vtkTextMapper
{
public:
  static vtkOpenGLTextMapper* New();
  vtkTypeMacro(vtkOpenGLTextMapper, vtkTextMapper);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  void RenderOverlay(vtkViewport* vp, vtkActor2D* act) override;

protected:
  vtkOpenGLTextMapper();
  ~vtkOpenGLTextMapper() override;

  void RenderGL2PS(vtkViewport* vp, vtkActor2D* act, vtkOpenGLGL2PSHelper* gl2ps);

private:
  vtkOpenGLTextMapper(const vtkOpenGLTextMapper&) = delete;
  void operator=(const vtkOpenGLTextMapper&) = delete;
};

#endif // vtkOpenGLTextMapper_h
