#!/usr/bin/env python
# coding: utf-8

import requests
from scipy import stats
import h5py
import numpy as np
from velocity import get

baseUrl = 'http://www.tng-project.org/api/'
headers = {"api-key":"47e1054245932c83855ab4b7af6a7df9"}


id = 25822
redshift = 2
scale_factor = 1.0 / (1+redshift)
little_h = 0.6774
solar_Z = 0.0127
url = "http://www.tng-project.org/api/TNG50-1/snapshots/z=" + str(redshift) + "/subhalos/" + str(id)

def particle_type(matter):

    if matter == "gas":
        part_type = 'PartType0'

    if matter == "stars":
        part_type = 'PartType4'

    params = {matter:'Coordinates,Masses'}

    sub = get(url)
    saved_filename = get(url+"/cutout.hdf5")

    with h5py.File(saved_filename,'r') as f:
        # NOTE! If the subhalo is near the edge of the box, you must take the 
        # periodic boundary into account! (we ignore it here)
        dx = f[part_type]['Coordinates'][:,0] - sub['pos_x']
        dy = f[part_type]['Coordinates'][:,1] - sub['pos_y']
        dz = f[part_type]['Coordinates'][:,2] - sub['pos_z']
        masses = f[part_type]['Masses'][:]*(10**10 / 0.6774)

        rr = np.sqrt(dx**2 + dy**2 + dz**2)
        rr *= scale_factor/little_h # ckpc/h -> physical kpc
        
        mass,bin_edge,num = stats.binned_statistic(rr,masses,statistic='sum',bins=np.linspace(0,30,50))
        
        f.close()
        
        return mass,bin_edge



def dm_mass():
    params = {'DM':'Coordinates,SubfindHsml'}

    sub = get(url)
    saved_filename = get(url+"/cutout.hdf5")

    with h5py.File(saved_filename,'r') as f:
        # NOTE! If the subhalo is near the edge of the box, you must take the 
        # periodic boundary into account! (we ignore it here)
        num = f['PartType1']['SubfindHsml'][:]

        dx = f['PartType1']['Coordinates'][:,0] - sub['pos_x']
        dy = f['PartType1']['Coordinates'][:,1] - sub['pos_y']
        dz = f['PartType1']['Coordinates'][:,2] - sub['pos_z']

        rr = np.sqrt(dx**2 + dy**2 + dz**2)
        rr *= scale_factor/little_h # ckpc/h -> physical kpc

        num_dm,bin_edge,x = stats.binned_statistic(rr,num,statistic='sum',bins=np.linspace(0,30,50))
        f.close()
        
        mass_dm_tot = 0.45*10**6
        mass_dm = num_dm*mass_dm_tot
        
        return mass_dm


    
    
def find_circ_vel():
    g = 'gas'
    stars = 'stars'

    mass_gas,r_gas = particle_type(g)
    mass_stars,r_stars = particle_type(stars)
    mass_dm = dm_mass() 
    
    r = (r_gas[1:]+r_gas[:-1])/2

    G_const = 4.30091*10**-6


    mass_enc_gas = np.cumsum(mass_gas)
    mass_enc_stars = np.cumsum(mass_stars)
    mass_enc_dm = np.cumsum(mass_dm)

    mass_tot = mass_enc_dm + mass_enc_gas + mass_enc_stars

    vel_circ = np.sqrt(G_const*mass_tot/r)
    
    return r,vel_circ


def star_pos_vel():
    params = {'stars':'Coordinates,Velocities'}

    sub = get(url)
    saved_filename = get(url+"/cutout.hdf5")
    
    with h5py.File(saved_filename,'r') as f:
        # NOTE! If the subhalo is near the edge of the box, you must take the 
        # periodic boundary into account! (we ignore it here)
        dx = (f['PartType4']['Coordinates'][:,0] - sub['pos_x'])*scale_factor
        dy = (f['PartType4']['Coordinates'][:,1] - sub['pos_y'])*scale_factor
        dz = (f['PartType4']['Coordinates'][:,2] - sub['pos_z'])*scale_factor

        vx = f['PartType4']['Velocities'][:,0]*np.sqrt(scale_factor) - sub['vel_x']
        vy = f['PartType4']['Velocities'][:,1]*np.sqrt(scale_factor) - sub['vel_y']
        vz = f['PartType4']['Velocities'][:,2]*np.sqrt(scale_factor) - sub['vel_z']

        star_masses = f['PartType4']['Masses'][:]*(10**10 / 0.6774)
        
        pos = np.array((dx,dy,dz)).T
        print(np.shape(pos))
        
        vel = np.array((vx,vy,vz)).T
        print(np.shape(vel))
        
    return(pos,vel,star_masses)