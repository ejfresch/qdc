#!/usr/bin/python

# Copyright 2009-2012 - Luca Freschi <l.freschi@gmail.com>                                                                       

# This file is part of QDC.                                                                                                  
# QDC is free software: you can redistribute it and/or modify                                                               
# it under the terms of the GNU General Public License as published by                                                      
# the Free Software Foundation, either version 3 of the License, or                                                         
# (at your option) any later version.  

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import sys
import os.path
import commands




def checks():
    msg=":: Preliminary checks..."
    if os.path.isdir('results'):
        status, out = commands.getstatusoutput('rm -rf results')
        if status !=0:
            msg=msg+"[ERROR]\n"
            return status, msg, out
    if os.path.isfile('engine'):
        status, out = commands.getstatusoutput('rm engine')
        if status !=0:
            msg=msg+"[ERROR]\n"
            return status, msg, out
    status, out = commands.getstatusoutput('mkdir results')
    if status !=0:
        msg=msg+"[ERROR]\n"
        return status, msg, out

    status=0
    out="[ok]\n"
    return status, msg, out


def compile_model(file):
    msg=":: Compilation of the model..."
    parse_file='./parser '+str(file)
    status, out = commands.getstatusoutput(parse_file)
    if status !=0:
        msg=msg+"[ERROR]\n"
        return status, msg, out
    status, out = commands.getstatusoutput('make engine')
    if status !=0:
        msg=msg+"[ERROR]\n"
        return status, msg, out
    status=0
    out="[ok]\n"
    return status, msg, out


def simulation():
    msg=":: Simulation..."
    status, out = commands.getstatusoutput('./engine '+' 0.1')
    if status !=0:
        msg=msg+"[ERROR]\n"
        return status, msg, out
    status=0
    out="[ok]\n"
    return status, msg, out



def write_results(i,model):
    msg=":: Output files..."
    
    base_name=os.path.basename(model)
    current=base_name+'_reagents'+str(i)+'.csv'
    cmd='mv '+model+'_reagents.csv '+'results/'+current
    status, out = commands.getstatusoutput(cmd)
    if status !=0:
        msg=msg+"[ERROR]\n"
        return status, msg, out
    current=base_name+'_reactions'+str(i)+'.csv'
    cmd='mv '+model+'_reactions.csv '+'results/'+current
    status, out = commands.getstatusoutput(cmd)
    if status !=0:
        msg=msg+"[ERROR]\n"
        return status, msg, out
    current=base_name+'_reactioncounts'+str(i)+'.csv'
    cmd='mv '+model+'_reactioncounts.csv '+'results/'+current
    status, out = commands.getstatusoutput(cmd)
    if status !=0:
        msg=msg+"[ERROR]\n"
        return status, msg, out
    current=base_name+'_log'+str(i)+'.txt'
    cmd='mv '+model+'_log.txt '+'results/'+current
    status, out = commands.getstatusoutput(cmd)
    if status !=0:
        msg=msg+"[ERROR]\n"
        return status, msg, out
    status=0
    out="[ok]\n"
    return status, msg, out


if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "usage: %s <file_input> <number_of_simulations>" % sys.argv[0]
        sys.exit()


    file=sys.argv[1]
    print file+"\n"
    #I cut the extension
    model=os.path.splitext(file)[0]
    print model+"\n"
    n_of_simulations=int(sys.argv[2])
    s, m, o =checks()
    print m+o
    if s !=0:
        sys.exit()
    s, m, o=compile_model(file)
    print m+o
    if s !=0:
        sys.exit()

    for i in range(n_of_simulations):
        print "Run "+str(i+1) 
        s, m, o=simulation()
        print m+o
        if s !=0:
            break
        s, m, o=write_results(i,model)
        print m+o
        if s!=0:
            break
   
        
