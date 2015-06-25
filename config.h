/**
 * @file  config.h
 * @brief config.h is the configuration file for c code in ReproBLAS.
 *
 * Feel free to modify anything here, but be sure to read the comments describing functionality
 */

/**
 * @brief The fold of double precision indexed types when fold is not specified.
 *
 * @author Peter Ahrens
 * @date   26 May 2015
 */
#define DIDEFAULTFOLD 3

/**
 * @brief The fold of single precision indexed types when fold is not specified.
 *
 * @author Peter Ahrens
 * @date   26 May 2015
 */
#define SIDEFAULTFOLD 3

/**
 * @brief The maximum double precision fold supported by the library.
 *
 * Increases to #DIMAXFOLD will allow users to use larger folds (meaning greater accuracy), but will also increase runtime (by a small amount) and the time it takes to autotune the library (by a large amount).
 *
 * @author Peter Ahrens
 * @date   26 May 2015
 */
#define DIMAXFOLD 4

/**
 * @brief The maximum single precision fold supported by the library.
 *
 * Increases to #SIMAXFOLD will allow users to use larger folds (meaning greater accuracy), but will also increase runtime (by a small amount) and the time it takes to autotune the library (by a large amount).
 *
 * @author Peter Ahrens
 * @date   26 May 2015
 */
#define SIMAXFOLD 4
