---
layout: default
title: Annular Pressure Modelling Tool
---

# Annular Pressure Modelling Tool

## Introduction

With the exception of the 'Background' and 'Assumptions' sections which will be updated as new critical background information and assumptions are identified, all included sections are presented in chronological order. 

## Background
In Horizontal Directional Drilling (HDD), drilling fluid is pumped downhole to cool and lubricate the drill bit, help extract cuttings, and create downhole stability. The gap between the borehole wall and the pipe installed within is referred to as the annulus. During drilling, the drilling fluid fills the annulus and exerts pressure on the ground surrounding the borehole, referred to as annular pressure. If the annular pressure becomes too high at a specific point, it can lead to uncontrolled soil fracture and subsequent escape of drilling fluid from the borehole, known as frac-out. Frac-out poses both environmental and infrastructure risks, and is considered a critical risk in HDD operation. To minimize the risk of frac-out, it is necessary to model annular pressure during operation.  

*IMAGE PLACEHOLDER: https://www.chasolutions.com/sites/default/assets/HDD_Blog/HDD-3a.jpg*

## Problem Statement


## Scope


## Assumptions

1. The bore path can be represented using a piecewise function of linear or arc-segment components

*IMAGE PLACEHOLDER: cred wikipedia user Krishnavedala*

2. The bore path is treated as 2-dimensional, with no left/right deviation

3. Borehole diameter *D<sub>b</sub>* and pipe outer diameter *D<sub>p</sub>* remain constant

4. The pipe is centered in the borehole (concentric annulus)

5. The drilling fluid is Newtonian and incompressible (constant dynamic viscosity *μ* and density *ρ*)

6. Steady flow

7. Flow is hydraulically smooth

8. The effect of drill cuttings on fluid properties is neglected 

## Core Functions
The first step is to begin creating the core functions required for the model. The first two functions (pictured below) calculate the cross-sectional area and hydraulic diameter of the concentric annulus, and are fairly self-explanatory.

*IMAGE PLACEHOLDER: corefunctions1*

Next, a path builder function is required to model the drill path based on user inputs. Due to the length of its code, no snapshots are included. The path builder approximates the drill path as a piecewise function of tangent and arc segments. Though 5 segments is standard, any number of segments can be used. The nature (tangent or arc) and vertical displacement of each segment are input by the user. Length and entry angle must be specified for tangent segments, while radius and deflection angle are needed for arc segments. 

*IMAGE PLACEHOLDER: create graphic for path modelling*



## References
