/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <DDqtControls.h>
#include <QAbstractSpinBox>
#include <QComboBox>
#include <QSlider>
#include <QAbstractScrollArea>
#include <QWheelEvent>
#include <QCoreApplication>

namespace model::controls
{

    bool WheelGuard::eventFilter(QObject* watched,QEvent* event)
    {
        if(event->type()!=QEvent::Wheel)
        {
            return false;
        }
        auto* widget=qobject_cast<QWidget*>(watched);
        QWidget* control=widget;
        while(control && !qobject_cast<QAbstractSpinBox*>(control)
              && !qobject_cast<QComboBox*>(control) && !qobject_cast<QSlider*>(control))
        {
            control=control->parentWidget();
        }
        if(!control)
        {
            return false;
        }
        auto* wheel=static_cast<QWheelEvent*>(event);
        for(auto* parent=control->parentWidget();parent;parent=parent->parentWidget())
        {
            if(auto* area=qobject_cast<QAbstractScrollArea*>(parent))
            {
                const QPointF local=area->viewport()->mapFromGlobal(wheel->globalPosition());
                QWheelEvent forwarded(local,wheel->globalPosition(),wheel->pixelDelta(),wheel->angleDelta(),
                                      wheel->buttons(),wheel->modifiers(),wheel->phase(),wheel->inverted(),wheel->source());
                QCoreApplication::sendEvent(area->viewport(),&forwarded);
                break;
            }
        }
        // Also suppress wheel edits when no surrounding scroll area exists.
        wheel->accept();
        return true;
    }

    auto GroupBox::showEvent(QShowEvent *event) -> void
    {
        if (isCheckable() && !body && layout())
        {
            body = new QWidget(this);
            body->setLayout(layout());
            auto *outer = new QVBoxLayout(this);
            outer->setContentsMargins(0, 0, 0, 0);
            outer->setSpacing(0);
            outer->addWidget(body);
            body->setVisible(isChecked());
        }
        QGroupBox::showEvent(event);
    }

    GroupBox::GroupBox(const QString &title, QWidget *parent) : QGroupBox(title, parent)
    {
        setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Maximum);
        connect(this, &QGroupBox::toggled, this,
                [this](bool checked)
                {
                    if (body)
                    {
                        body->setVisible(checked);
                        updateGeometry();
                    }
                });
    }
} // namespace model::controls
