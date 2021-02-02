/*
 * defines.h
 *
 *  Created on: 10 Feb 2020
 *      Author: jkiesele
 */

#ifndef MINICALO_INCLUDE_DEFINES_H_
#define MINICALO_INCLUDE_DEFINES_H_

//#define NOPYTHIA

//#define USEVIS

//#define USEPYTHIA B4PartGeneratorBase::minbias

#ifdef USEPYTHIA
#define ONLY_ENERGY_OUTPUT
#endif

//#define FROMBEAMSPOT

#endif /* MINICALO_INCLUDE_DEFINES_H_ */
