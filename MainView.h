#ifndef CITYVIEW_H
#define CITYVIEW_H

#include <QGraphicsView>

class MainView : public QGraphicsView
{
    Q_OBJECT

public:
    MainView(QWidget *parent = 0);

protected:
	void wheelEvent(QWheelEvent *event);
	// can not drag the image if using this function
// 	void mousePressEvent(QMouseEvent * event);

signals :	
	void mousePressed(int x,int y);

};

#endif
