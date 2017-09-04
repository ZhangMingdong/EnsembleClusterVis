#include <QtGui>
#include <cmath>

#include <QDebug>
#include <QPoint>

#include "MainView.h"

MainView::MainView(QWidget *parent)
    : QGraphicsView(parent)
{
    setDragMode(ScrollHandDrag);
}

void MainView::wheelEvent(QWheelEvent *event)
{
    double numDegrees = -event->delta() / 8.0;
    double numSteps = numDegrees / 15.0;
    double factor = std::pow(1.125, numSteps);
    scale(factor, factor);
}

// void MainView::mousePressEvent(QMouseEvent * event){
// 	QPoint pt = event->pos();
// 	emit mousePressed(pt.x(),pt.y());
// }