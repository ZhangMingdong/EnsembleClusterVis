#pragma once

#include <QWidget>
#include "ui_ControlWidget.h"

class ControlWidget : public QWidget
{
	Q_OBJECT

public:
	ControlWidget(QWidget *parent = Q_NULLPTR);
	~ControlWidget();

public:
	Ui::ControlWidget ui;
};
