#pragma once

#ifndef SIMUTIL_DIMENSION_H
#define SIMUTIL_DIMENSION_H


#ifdef THREE_DIMENSION
#define DIMENSION 3
#define Vec Vec3
#define Vecd Vec3d
#else
#define DIMENSION 2
#define Vec Vec2
#define Vecd Vec2d
#endif // THREE_DIMENSION


#endif // !DIMENSION_H
