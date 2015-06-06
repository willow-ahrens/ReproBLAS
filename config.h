/**
 * @file  config.h
 * @brief config.h is the configuration file for c code in ReproBLAS.
 *
 * Feel free to modify anything here, but be sure to read the comments describing functionality
 */

/**
 * @brief The fold of indexed types used anywhere fold is not specified.
 *
 * @author Peter Ahrens
 * @date   26 May 2015
 */
#define DEFAULT_FOLD 3

/**
 * @brief The maximum fold supported by the library.
 *
 * Increases to #MAX_FOLD will allow users to use larger folds (meaning greater accuracy), but will also increase suntime (by a small amount) and the time it takes to autotune the library (by a large amount).
 *
 * @author Peter Ahrens
 * @date   26 May 2015
 */
#define MAX_FOLD 4

/**
 * @brief The maximum number of carry overflow supported by the library.
 *
 * Increases to #MAX_FLOW will allow users to use more carries per fold.
 * Future support for larger values of MAX_FLOW is being considered, but
 * currently the only supported value is 1.
 *
 * @author Peter Ahrens
 * @date   5 Jun 2015
 */
#define MAX_FLOW 1
