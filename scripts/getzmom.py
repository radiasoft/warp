from warp import *
getzmom_version = "$Id: getzmom.py,v 1.21 2006/04/28 22:20:53 dave Exp $"

def getzmomdoc():
  print """
zmmnt  makes appropriate calls to compiled code to calculate the
       particle moments
  """

def zmmnt(itask=0,js=None,jslist=None,groupsize=256):
  """
zmmnt(itask=0,js=None, jslist=xrange(0,top.ns))
  makes appropriate calls to compiled code to calculate the
  particle moments
  - itask=0 task to carry out
    0 All steps to calculate moments (doing steps 1, 2, and 3 below)
    1 initialization, zeros out moments variables
    2 gather moments
    3 final processing of moments, averaging and derived moments
  - js: select species over which to calculate moments,
        None: calculates moments over range 0 .. top.ns or,
        alternatively, over list of moments specified in jslist
  - jslist: list of species to calculate moments for
  - groupsize=256: size of particle groups sent to getzmmnt
                   Use caution since temp variables are used that are
                   18*groupsize in size.
  """

  # If changed, fix Win_Moments arrays
  if (len(top.pnumz) != top.nzmmnt+1): gchange ("Z_Moments")

  ismax = maxval(top.pgroup.sid)+1
  
  # --- This is bad idea and so is commented out.
  # If itask is greater than zero, assume that the intent was to do the
  # calculation regardless of the value of laccumulate_zmoments.  Save
  # the value of laccumulate_zmoments and then set it to false.
# if (itask > 0):
#   l1 = top.laccumulate_zmoments
#   top.laccumulate_zmoments = false

  # Zero out the moments arrays
  if (itask == 0 or itask == 1):
    getzmmnt(0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,1,0,0.,0.,0.,1,1,top.ns,
             top.tempmaxp,top.tempminp,top.tempzmmnts0,top.tempzmmnts)

  # Calculate the moments
  if jslist is None:
    if js is None:
      jslist = xrange(top.pgroup.ns)
    else:
      jslist = [js]
  if (itask == 0 or itask == 2):
    for js in jslist:
      isid = top.pgroup.sid(js) + 1
      for ipmin in xrange(top.pgroup.ins[js]-1,top.pgroup.ins[js]+top.pgroup.nps[js]-1,groupsize):
         ip = min(groupsize, top.pgroup.ins[js]+top.pgroup.nps[js]-ipmin-1)
         try:
           weighted = top.wpid
         except AttributeError:
           weighted = 0
         x = top.pgroup.xp[ipmin:ipmin+ip]
         y = top.pgroup.yp[ipmin:ipmin+ip]
         z = top.pgroup.zp[ipmin:ipmin+ip]
         ux = top.pgroup.uxp[ipmin:ipmin+ip]
         uy = top.pgroup.uyp[ipmin:ipmin+ip]
         uz = top.pgroup.uzp[ipmin:ipmin+ip]
         gaminv = top.pgroup.gaminv[ipmin:ipmin+ip]
         if(not weighted):
           getzmmnt(ip,x,y,z,ux,uy,uz,gaminv,
                    top.pgroup.sq[js],top.pgroup.sm[js],top.pgroup.sw[js],
                    top.dt,top.pgroup.dtscale[js],2,top.nplive,ux,uy,uz,
                    js+1,isid,ismax,
                    top.tempmaxp,top.tempminp,top.tempzmmnts0,top.tempzmmnts)
         else:
           pid = top.pgroup.pid[ipmin:ipmin+ip,top.wpid-1]
           getzmmnt_weights(ip,x,y,z,ux,uy,uz,gaminv,pid,
                    top.pgroup.sq[js],top.pgroup.sm[js],top.pgroup.sw[js],
                    top.dt,top.pgroup.dtscale[js],2,top.nplive,
                    ux,uy,uz,js+1,isid,ismax,
                    top.tempmaxp,top.tempminp,top.tempzmmnts0,top.tempzmmnts)

  # Do final summing and averaging of the moments
  if (itask == 0 or itask == 3):
    getzmmnt(0,0.,0.,0.,0.,0.,0.,0.,
             0.,0.,0.,0.,0.,3,top.nplive,0.,0.,0.,1,1,top.ns,
             top.tempmaxp,top.tempminp,top.tempzmmnts0,top.tempzmmnts)

  # Restore the value of laccumulate_zmoments
# if (itask > 0):
#   top.laccumulate_zmoments = l1

