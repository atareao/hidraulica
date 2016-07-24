#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# This file is part of Hidraulica
#
# Copyright (C) 2016 Lorenzo Carbonell
# lorenzo.carbonell.cerezo@gmail.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import uno
import math
import numpy
import unohelper
from es.atareao.libreoffice.Hidraulica import XHidraulica

# Template
# Template created by Lorenzo Carbonell <lorenzo.carbonell.cerezo@gmail.com>


class HidraulicaImpl(unohelper.Base, XHidraulica):
    def __init__(self, ctx):
        self.ctx = ctx

# -------------------- Functions to modify -----------------------------

    def reynolds(self, densidad, velocidad, diametro, viscosidaddinamica):
        ans = densidad * velocidad * diametro / viscosidaddinamica
        return ans

    def reynoldsaguadulce(self, temperatura, velocidad, diametro):
        viscosidaddinamicaagua = self.viscosidaddinamicaagua(temperatura)
        densidadagua = self.densidadagua(temperatura, 0.5, 1)
        return self.reynolds(densidadagua, velocidad, diametro,
                             viscosidaddinamicaagua)

    def reynoldsaguademar(self, temperatura, velocidad, diametro):
        viscosidaddinamicaagua = self.viscosidaddinamicaagua(temperatura)
        densidadagua = self.densidadagua(temperatura, 40, 1)
        return self.reynolds(densidadagua, velocidad, diametro,
                             viscosidaddinamicaagua)

    def viscosidaddinamicaagua(self, t):
        ans = (1.78757149000000E-03 * math.pow(t, 0) -
               5.86440000000000E-05 * math.pow(t, 1) +
               1.26529900857808E-06 * math.pow(t, 2) -
               1.70874544180172E-08 * math.pow(t, 3) +
               1.25746405762044E-10 * math.pow(t, 4) -
               3.78265790976119E-13 * math.pow(t, 5))
        return ans

    def densidadagua(self, T, s, p):
        # https://es.wikipedia.org/wiki/\
        # Anexo:C%C3%A1lculo_de_la_densidad_del_agua_del_mar
        # Función Densidad del océano, la cual la calcula a partir de T(°C),
        # s(psu) y p(bar)
        # La funcion rho(T,s,p) calcula la densidad del agua de mar
        # a partir de la aproximacion empirica de UNESCO del año 1981
        # Utilice T(Celsius), s(psu), p(bar)
        # Salida en unidades SI [kg/m^3]
        base = numpy.array([math.pow(T, 0), math.pow(T, 1), math.pow(T, 2),
                            math.pow(T, 3), math.pow(T, 4), math.pow(T, 5)])

        A = numpy.dot(numpy.array(
            [999.8425, 6.7939e-2, -9.0952e-3, 1.0016e-4, -1.12e-6, 6.53e-9]),
            base)
        B = numpy.dot(numpy.array(
            [8.2449e-1, -4.0899e-3, 7.6438e-5, -8.2467e-7, 5.3875e-9, 0]),
            base)
        C = numpy.dot(numpy.array(
            [-5.7246e-3, 1.0227e-4, -1.6546e-6, 0, 0, 0]),
            base)
        D = numpy.dot(numpy.array(
            [4.8314e-4, 0, 0, 0, 0, 0]),
            base)
        if p == 0:
            rho = A + B*s + C*math.pow(s, 1.5) + D*math.pow(s, 2)
        else:
            kt = self.Kt(T, s, p)
            if kt > 0:
                rho = (A + B*s + C*math.pow(s, 1.5) + D*math.pow(s, 2))/(
                       1-(p / kt))
        return rho

    def Kt(self, T, s, p):
        # Función Modulo de Compresibilidad Secante
        # Calcula el polinomio usando los parámetros entregados, y los envía a
        # la función rho(T,s,p)
        base = numpy.array([math.pow(T, 0), math.pow(T, 1), math.pow(T, 2),
                            math.pow(T, 3), math.pow(T, 4), math.pow(T, 5)])
        E = numpy.dot(numpy.array(
            [19652.21, 148.4206, -2.3271, 1.3604e-2, -5.1552e-5, 0]), base)
        F = numpy.dot(numpy.array(
            [54.6746, -0.6034, 1.0998e-2, -6.1670e-5, 0, 0]), base)
        G = numpy.dot(numpy.array(
            [7.944e-2, 1.6483e-2, -5.3009e-4, 0, 0, 0]), base)
        H = numpy.dot(numpy.array(
            [3.2399, 1.4371e-3, 1.1609e-4, -5.7790e-7, 0, 0]), base)
        I = numpy.dot(numpy.array(
            [2.2838e-3, -1.0981e-5, -1.6078e-6, 0, 0, 0]), base)
        J = numpy.dot(numpy.array(
            [1.9107e-4, 0, 0, 0, 0, 0]), base)
        M = numpy.dot(numpy.array(
            [8.5093e-5, -6.1229e-6, 5.2787e-7, 0, 0, 0]), base)
        N = numpy.dot(numpy.array(
            [-9.9348e-7, 2.0816e-8, 9.1697e-10, 0, 0, 0]), base)
        ans = (E + F*s + G*math.pow(s, 1.5) + (H + I*s +
               J*math.pow(s, 1.5))*p + (M + N*s)*math.pow(p, 2))
        return ans

    def friccion(self, diametro, caudal, rugosidad, reynolds):
        # diametro -> m
        # caudal -> m3/s
        # rugosidad -> mm
        # inicializaciones
        contador = 0
        f0 = 0.1
        f1 = 0.01
        # calculos
        velocidad = 4.0 * caudal / (math.pi * math.pow(diametro, 2.0))
        while abs((f1 - f0) / f0) > 0.00001 and contador < 100:
            contador = contador + 1
            f0 = f1
            f1 = 1 / math.pow(
                2 * math.log10(
                    rugosidad / (3.71 * diametro * 1000.0) + 2.51 / (
                        reynolds * math.sqrt(f0))), 2.0)
        return f0

    def friccionaguademar(self, diametro, caudal, rugosidad,
                             temperatura):
        velocidad = 4.0 * caudal / (math.pi * math.pow(diametro, 2.0))
        reynolds = self.reynoldsaguademar(temperatura, velocidad,
                                             diametro)
        return self.friccion(diametro, caudal, rugosidad, reynolds)

    def friccionaguadulce(self, diametro, caudal, rugosidad,
                            temperatura):
        velocidad = 4.0 * caudal / (math.pi * math.pow(diametro, 2.0))
        reynolds = self.reynoldsaguadulce(temperatura, velocidad,
                                            diametro)
        return self.friccion(diametro, caudal, rugosidad, reynolds)

    def perdidacarga(self, friccion, longitud, caudal, diametro):
        return 8.0 * friccion * longitud * math.pow(caudal, 2.0) / (
            math.pow(math.pi, 2.0) * 9.81 * math.pow(diametro, 5.0))

    def velocidad(self, caudal, diametro):
        return 4.0 * caudal / (math.pi * math.pow(diametro, 2.0))


def createInstance(ctx):
    return HidraulicaImpl(ctx)

g_ImplementationHelper = unohelper.ImplementationHelper()
g_ImplementationHelper.addImplementation(
    createInstance, "es.atareao.libreoffice.Hidraulica.python.HidraulicaImpl",
    ("com.sun.star.sheet.AddIn",),)

if __name__ == '__main__':
    fi = HidraulicaImpl(None)
    print(fi.sample())
    print(fi.mult25(3, 2))
    print(fi.densidadagua(94, 0.5, 1))
    print(fi.densidadagua(0, 0, 1))
    print(fi.densidadagua(20, 0, 1))
    print(fi.densidadagua(20, 0, 2))
    print(fi.densidadagua(20, 0, 3))
    print(fi.densidadagua(20, 0, 4))
    print(fi.densidadagua(20, 0, 5))
    print(fi.densidadagua(20, 0, 6))
    print(fi.reynoldsaguadulce(20, 1, 1))
    print(fi.reynoldsaguademar(20, 1, 1))
    diametro = 1
    caudal = 2
    rugosidad = 3
    temperatura = 20
    velocidad = 4.0 * caudal / (math.pi * math.pow(diametro, 2.0))
    print(velocidad)
    print(fi.friccionaguadulce(diametro, caudal, rugosidad, temperatura))
    print(fi.friccionaguademar(diametro, caudal, rugosidad, temperatura))
    print(fi.perdidacarga(0.02, 1000, caudal, diametro))
