/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Charts module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#include "customtablemodel.h"
#include "def.h"


#include <QtCore/QVector>
#include <QtCore/QTime>
#include <QtCore/QRect>
#include <QtGui/QColor>

#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QDebug>


void CustomTableModel::SetIndex(int nIndex){
// 	m_data.clear();
// 	for (int i = 0; i < m_rowCount; i++) {
// 		QModelIndex index(i,0);
// 		index.row()
// 
// 		QVector<qreal>* dataVec = new QVector<qreal>(m_columnCount);
// 		dataVec->replace(0, i / 10.0 - 50);
// 		dataVec->replace(1, _pData[nIndex*1001+i]);
// 		m_data.append(dataVec);
// 	}
}
CustomTableModel::CustomTableModel(double* pData, QObject *parent) :
    QAbstractTableModel(parent)
{
	qsrand(QTime(0, 0, 0).secsTo(QTime::currentTime()));

	_pData = pData;
	m_data.clear();

	m_columnCount = 2;
	m_rowCount = 1401;

	// m_data
	for (int i = 0; i < m_rowCount; i++) {
		QVector<qreal>* dataVec = new QVector<qreal>(m_columnCount);
		dataVec->replace(0, i/10.0 - 70);
		dataVec->replace(1, _pData[i]);

// 		for (int j = 0; j < m_columnCount - 1; j++)
// 		{
// 			dataVec->replace(j + 1, _pData[j*1001+i]);
// 		}
		m_data.append(dataVec);
	}

}

int CustomTableModel::rowCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent)
    return m_data.count();
}

int CustomTableModel::columnCount(const QModelIndex &parent) const
{
    Q_UNUSED(parent)
    return m_columnCount;
}

QVariant CustomTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (role != Qt::DisplayRole)
        return QVariant();

    if (orientation == Qt::Horizontal) {
        if (section % 2 == 0)
            return "x";
        else
            return "y";
    } else {
        return QString("%1").arg(section + 1);
    }
}

QVariant CustomTableModel::data(const QModelIndex &index, int role) const
{
    if (role == Qt::DisplayRole) {
        return m_data[index.row()]->at(index.column());
    } else if (role == Qt::EditRole) {
        return m_data[index.row()]->at(index.column());
    } else if (role == Qt::BackgroundRole) {
        foreach (QRect rect, m_mapping) {
            if (rect.contains(index.column(), index.row()))
                return QColor(m_mapping.key(rect));
        }
        // cell not mapped return white color
        return QColor(Qt::white);
    }
    return QVariant();
}

bool CustomTableModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
    if (index.isValid() && role == Qt::EditRole) {
        m_data[index.row()]->replace(index.column(), value.toDouble());
        emit dataChanged(index, index);
        return true;
    }
    return false;
}

Qt::ItemFlags CustomTableModel::flags(const QModelIndex &index) const
{
    return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
}

void CustomTableModel::addMapping(QString color, QRect area)
{
    m_mapping.insertMulti(color, area);
}
