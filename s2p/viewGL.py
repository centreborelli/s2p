#!/usr/bin/env python
# Copyright (C) 2015, Gabriele Facciolo <gfacciol@gmail.com>

## {{{ http://code.activestate.com/recipes/325391/ (r1)
from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
import sys
from PIL import Image

name = 'select a region'

# global variables for the mouse
x0=0; y0=0; w0=0; h0=0; b0state='';

# global variables for the image size
ix,iy=100,100
# and for the relative displacement
dx,dy=0,0
imageData=0


def loadImage(imageName):
    im = Image.open(imageName)
    try:
       ix, iy, image = im.size[0], im.size[1], im.convert("RGBA").tobytes("raw", "RGBA", 0, -1)
    except SystemError:
       ix, iy, image = im.size[0], im.size[1], im.convert("RGBA").tobytes("raw", "RGBX", 0, -1)
    return image, ix, iy



def display( ):
    """Render scene geometry"""

    glClear(GL_COLOR_BUFFER_BIT);

#    glDisable( GL_LIGHTING) # context lights by default

    # this is the effective size of the current window
    # ideally winx,winy should be equal to ix,iy BUT if the
    # image is larger than the screen glutReshapeWindow(ix,iy)
    # will fail and winx,winy will be truncated to the size of the screen
    winx=glutGet(GLUT_WINDOW_WIDTH)
    winy=glutGet(GLUT_WINDOW_HEIGHT)


    # setup camera
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho (0, winx, winy, 0, -1, 1);


    glEnable (GL_TEXTURE_2D); #/* enable texture mapping */
    glBindTexture (GL_TEXTURE_2D, 13); #/* bind to our texture, has id of 13 */

    global dx,dy,ix,iy

    glBegin( GL_QUADS );
    glColor3f(1.0, 0.0, 0.0);
    glTexCoord2d(0.0,0.0); glVertex3d(0 -dx,iy-dy,0);
    glTexCoord2d(1.0,0.0); glVertex3d(ix-dx,iy-dy,0);
    glTexCoord2d(1.0,1.0); glVertex3d(ix-dx,0 -dy,0);
    glTexCoord2d(0.0,1.0); glVertex3d(0 -dx,0 -dy,0);
    glEnd();
    glDisable (GL_TEXTURE_2D); #/* disable texture mapping */


    # show region
    global x0,y0,w0,h0,b0state
    if(b0state=='pressed'):
       glEnable (GL_BLEND)
       glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

       glBegin( GL_QUADS );
       glColor4f(1.0, 0.0, 0.0,.6);
       glVertex3d(x0   ,y0   ,-0.1);
       glColor4f(0.0, 1.0, 1.0,.3);
       glVertex3d(x0+w0,y0   ,-0.1);
       glColor4f(0.0, 0.0, 1.0,.6);
       glVertex3d(x0+w0,y0+h0,-0.1);
       glColor4f(0.0, 1.0, 0.0,.3);
       glVertex3d(x0   ,y0+h0,-0.1);
       glEnd();

       glDisable (GL_BLEND)


    glFlush()
    glutSwapBuffers()


def setupTexture(imageData, ix,iy):
    """texture environment setup"""
    glBindTexture(GL_TEXTURE_2D, 13)
    glPixelStorei(GL_UNPACK_ALIGNMENT,1)
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
#    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
#    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL)

    glTexImage2D(
      GL_TEXTURE_2D, 0, GL_RGBA, ix, iy, 0,
      GL_RGBA, GL_UNSIGNED_BYTE, imageData
      )


def mouseMotionPass(x,y):
    global dx,dy
    title='p:%s,%s' % (x+dx,y+dy)
    glutSetWindowTitle(title)


def mouseMotion(x,y):
    global x0,y0,w0,h0,b0state
    global dx,dy
    w0,h0 = x-x0,y-y0
    title='p:%s,%s [+%s+%s %sx%s]' % (x+dx,y+dy,x0+dx,y0+dy,w0,h0)
    glutSetWindowTitle(title)
    glutPostRedisplay()


def mouseButtons(button, state, x,y):
    global x0,y0,w0,h0,b0state
    if button==GLUT_LEFT_BUTTON and state==GLUT_DOWN:
       x0,y0=x,y
       b0state='pressed'
    elif button==GLUT_LEFT_BUTTON and state==GLUT_UP:
       w0,h0 = x-x0,y-y0
       b0state='released'
       print x0+dx, y0+dy, x+dx, y+dy
       exit(0)



# letters and numbers
def keyboard(key, x, y):
    if key=='q':
       exit(0)
#    print ord(key)



# handle arrow keys
def keyboard2(key, x, y):
    global dx,dy
    winx=glutGet(GLUT_WINDOW_WIDTH)
    winy=glutGet(GLUT_WINDOW_HEIGHT)
    if key==102:
       dx=dx+int(winx/16)
    elif key==101:
       dy=dy-int(winy/16)
    elif key==100:
       dx=dx-int(winx/16)
    elif key==103:
       dy=dy+int(winy/16)
#    print key
    glutPostRedisplay()




def main():
    # verify input
    if len(sys.argv) > 1:
       I1 = sys.argv[1]
    else:
       print "Incorrect syntax, use:"
       print '  > ' + sys.argv[0] + " image.png"
       # show default image if exists
       I1 = '/Users/facciolo/uiskentuie_standing_stone.png'
       try:
          os.stat(I1)
       except OSError:
          exit(1)

    # globals
    global ix,iy,imageData

    #init opengl
    glutInit(sys.argv)
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH)
    glutInitWindowSize(ix,iy)
    glutCreateWindow(name)
    glutDisplayFunc(display)
#    glutReshapeFunc (glutResize);
#    glutIdleFunc(display)
    glutMouseFunc(mouseButtons)
    glutKeyboardFunc(keyboard)
    glutSpecialFunc(keyboard2)
    glutPassiveMotionFunc(mouseMotionPass);
    glutMotionFunc(mouseMotion);

    # read the image
    (imageData,ix,iy) = loadImage(I1)

    glutReshapeWindow(ix,iy)

    # setup texture
    setupTexture(imageData,ix,iy)

    glutMainLoop()


if __name__ == '__main__': main()
## end of http://code.activestate.com/recipes/325391/ }}}
