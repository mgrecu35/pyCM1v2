
def write_sub(var_def):
    i1=var_def.find("(")
    var_name=var_def[:i1]
    var_dims=var_def[i1:]
    print(var_name,var_dims)
    l1='subroutine set_%s(ibmo,iemo,jbmo,jemo,kbmo,kemo,%s_in)\n'%(var_name,\
                                                                 var_name)
    l2='use cm1vars\n'
    l3='implicit none\n'
    l4='integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo\n'
    l5='real, intent(in) :: %s%s\n'%(var_name+'_in',var_dims)
    for dim1 in ["ib","ie","jb","je","kb","ke"]:
        l5=l5.replace(dim1,dim1+"mo")
    l6='%s=%s_in\n'%(var_name,var_name)
    l7='end subroutine set_%s\n'%var_name
    return l1+l2+l3+l4+l5+l6+l7

def write_get_sub(var_def):
    i1=var_def.find("(")
    var_name=var_def[:i1]
    var_dims=var_def[i1:]
    print(var_name,var_dims)
    l1='subroutine get_%s(ibmo,iemo,jbmo,jemo,kbmo,kemo,%s_out)\n'%(var_name,\
                                                                 var_name)
    l2='use cm1vars\n'
    l3='implicit none\n'
    l4='integer :: ibmo,iemo,jbmo,jemo,kbmo,kemo\n'
    l5='real, intent(out) :: %s%s\n'%(var_name+'_out',var_dims)
    for dim1 in ["ib","ie","jb","je","kb","ke"]:
        l5=l5.replace(dim1,dim1+"mo")
    l6='%s_out=%s\n'%(var_name,var_name)
    l7='end subroutine get_%s\n'%var_name
    return l1+l2+l3+l4+l5+l6+l7


var_def='u3d(ib:ie+1,jb:je,kb:ke)'

code=write_sub(var_def)

vars_def=['u3d(ib:ie+1,jb:je,kb:ke)','v3d(ib:ie,jb:je+1,kb:ke)',\
          'w3d(ib:ie,jb:je,kb:ke+1)','ua(ib:ie+1,jb:je,kb:ke)',\
          'va(ib:ie,jb:je+1,kb:ke)','wa(ib:ie,jb:je,kb:ke+1)',\
          'ppi(ib:ie,jb:je,kb:ke)','pp3d(ib:ie,jb:je,kb:ke)']

#'tha(ib:ie,jb:je,kb:ke)','th3d(ib:ie,jb:je,kb:ke)'
from define_set import *
f=open("../src/cm1var_setters.f90","w")
for var_def in vars_def:
    code=write_sub(var_def)
    for l in code:
        f.write(l)
    f.write("!------------------\n")
f.write("!---------Getters--------\n")

for var_def in vars_def:
    code=write_get_sub(var_def)
    for l in code:
        f.write(l)
    f.write("!------------------\n")

f.close()

