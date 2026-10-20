/* This file is part of MoDELib, distributed under GNU GPL v2. */
#include <DDqtControls.h>

namespace model::controls
{

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
