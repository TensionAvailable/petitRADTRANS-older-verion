C     -*- fortran -*-
C     This file is autogenerated with f2py (version:1.20.3)
C     It contains Fortran 77 wrappers to fortran functions.

      subroutine f2pywrapinteg_parab (integ_parabf2pywrap, x, y, z
     &, fx, fy, fz, a, b)
      external integ_parab
      double precision x
      double precision y
      double precision z
      double precision fx
      double precision fy
      double precision fz
      double precision a
      double precision b
      double precision integ_parabf2pywrap, integ_parab
      integ_parabf2pywrap = integ_parab(x, y, z, fx, fy, fz, a, b)
      end


      subroutine f2pywraphansen_size_nr (hansen_size_nrf2pywrap, r
     &, a, b)
      external hansen_size_nr
      double precision r
      double precision a
      double precision b
      double precision hansen_size_nrf2pywrap, hansen_size_nr
      hansen_size_nrf2pywrap = hansen_size_nr(r, a, b)
      end


      subroutine f2pywraphansen_size_dndr (hansen_size_dndrf2pywra
     &p, r, a, b, k)
      external hansen_size_dndr
      double precision r
      double precision a
      double precision b
      double precision k
      double precision hansen_size_dndrf2pywrap, hansen_size_dndr
      hansen_size_dndrf2pywrap = hansen_size_dndr(r, a, b, k)
      end


      subroutine f2pywrapbisect_particle_rad (bisect_particle_radf
     &2pywrap, x1, x2, gravity, rho, rho_p, temp, mmw, w_star)
      external bisect_particle_rad
      double precision x1
      double precision x2
      double precision gravity
      double precision rho
      double precision rho_p
      double precision temp
      double precision mmw
      double precision w_star
      double precision bisect_particle_radf2pywrap, bisect_particl
     &e_rad
      bisect_particle_radf2pywrap = bisect_particle_rad(x1, x2, gr
     &avity, rho, rho_p, temp, mmw, w_star)
      end

