"""Graphical tuning of lattice elements with the envelope code"""

from warp import *
envtuner_version = "$Id: envtuner.py,v 1.3 2009/09/14 23:27:17 dave Exp $"

def makeenvtuneplot(zl,zu,qz):
    # --- Make the envelope plot, setting and returning the plot limits.
    plg(env.aenv,env.zenv)
    plg(env.benv,env.zenv,color='red')
    #plg(env.apenv,env.zenv,color='red')
    envamin = min(env.aenv)
    envbmin = min(env.benv)
    envmin = min(envamin,envbmin)
    envamax = max(env.aenv)
    envbmax = max(env.benv)
    envmax = max(envamax,envbmax)
    limits() # --- reset the limits to the extremum of the data
    p = limits() # --- get the plot limits
    # --- Plot vertical guiding lines between elements
    for z in qz:
        plg(array([p[2],p[3]]),array([z,z]),type='dot')
    return p


def envtuner(elem='quad',suffix='de',zl=env.zl,zu=env.zu,
             scale=1.,mindelta=1.e-3):
    """In the gist plotting window, the mouse is used to modify
the lattice elements.
elem='quad': element type to modify
suffix='de': suffix of quantity to modify
zl=env.zl: start of range to modify
zu=env.zu: end of range to modify
scale=1.: scale factor to convert the mouse movement to the field strength.
mindelta=1.e-3: minimum change to use if there was no mouse movement during
                the mouse click. If the mouse click is made in the lower
                portion of the plot, then apply -mindelta, otherwise +mindelta.

Left mouse button is used to scale the element. Click and drag between the
vertical lines to scale the element that is between the lines. The quantity is
multiplied by (1+scale*delta), where scale is the input argument, and delta is
the vertical fraction of the plot that the mouse was dragged over. An upward
drag gives a positive delta and downward a negative delta. Clicking without
dragging the mouse gives a value of mindelta, with a positive value if clicked
in the upper half of the plot and negative in the lower half.

Left mouse button with control key pressed changes the sign of the quantity.

Middle mouse button clears the frame.
Middle mouse button with control resets the quantity to its initial value.

Right mouse button exits.

    """
    nelem = getattr(top,'n'+elem)
    elemzs = getattr(top,elem+'zs')
    elemze = getattr(top,elem+'ze')
    elemqu = getattr(top,elem+suffix)

    # --- Find the lattice elements within the range.
    ielems = 0
    ieleme = 0
    while ielems < nelem and elemze[ielems]+top.zlatstrt < zl:
        ielems = ielems + 1
    while ieleme < nelem and elemzs[ieleme]+top.zlatstrt < zu:
        ieleme = ieleme + 1
    qz = zeros(ieleme - ielems + 1,'d')
    qz[1:] = (top.zlatstrt +
              0.5*(elemze[ielems:ieleme] + elemzs[ielems+1:ieleme+1]))
    if qz[-1] > zu: qz[-1] = zu
    if ielems == 0:
        qz[0] = zl
    else:
        qz[0] = top.zlatstrt + 0.5*(elemze[ielems-1] + elemzs[ielems])

    # --- Save the initial value so that the quantity can be reset.
    elemquinit = elemqu[ielems:ieleme+1].copy()

    # --- Make initial plot
    p = makeenvtuneplot(zl,zu,qz)

    # --- Now, enter the loop.
    done = 0
    while not done:

        # --- Get the mouse motion
        m = mouse(1,2,'')

        # --- If null was returned or if right button was clicked, the quit.
        if not m or m[9] == 3:
            done = 1
            continue

        # --- If middle button if pushed, then either clear the plot
        # --- or reset the data to the initial values.
        elif m[9] == 2:
            if m[10] == 0:
                fma()
                p = makeenvtuneplot(zl,zu,qz)
                continue
            elif m[10] == 4:
                elemqu[ielems:ieleme+1] = elemquinit
                delta = 0.

        elif m[9] == 1:
            # --- Get z position and which element it is in.
            z = m[0]
            for iq in range(len(qz)-1):
                if qz[iq] < z and z < qz[iq+1]: break
            if z > qz[iq+1]: iq += 1

            if m[10] == 0:
                # --- Get value of delta requested.
                delta = (m[3] - m[1])/(p[3] - p[2])
                if delta == 0:
                    if m[1] < 0.5*(p[2]+p[3]):
                        delta = -mindelta
                    else:
                        delta = +mindelta
            elif m[10] == 4:
                # --- With the control key, change the sign of the quantity
                delta = -2./scale

        # --- Change element quantity appropriately and print out the new value.
        elemqu[ielems+iq] *= (1. + delta*scale)
        print elem+suffix+' = '+repr(elemqu[ielems:ieleme+1])

        # --- Recalculate envelope and redo the plot.
        step()
        p = makeenvtuneplot(zl,zu,qz)
