/* This file is part of MoDELib, distributed under GNU GPL v2. */
#ifndef model_DDqtControls_H_
#define model_DDqtControls_H_
#include <QCheckBox>
#include <QGroupBox>
#include <QGridLayout>
#include <QLabel>
#include <QVBoxLayout>
#include <QShowEvent>

namespace model::controls
{
    // Keep the contents in one widget so collapsing preserves each child's own
    // visibility, enabled state, and nested group selection.
    class GroupBox : public QGroupBox
    {
        QWidget *body = nullptr;

      protected:
        void showEvent(QShowEvent *event) override
        {
            if(isCheckable() && !body && layout())
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

      public:
        explicit GroupBox(const QString &title, QWidget *parent = nullptr) : QGroupBox(title, parent)
        {
            setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Maximum);
            connect(this, &QGroupBox::toggled, this, [this](bool checked)
            {
                if(body)
                {
                    body->setVisible(checked);
                    updateGeometry();
                }
            });
        }
    };

    inline QGroupBox *box(const QString &title, QWidget *parent, const QString &description,
                          QCheckBox *toggle = nullptr)
    {
        auto *group = new GroupBox(title, parent);
        group->setToolTip(description);
        group->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Maximum);
        new QGridLayout(group);
        if(toggle)
        {
            // Preserve existing actor handlers and programmatic state changes;
            // the visible group checkbox is the single user-facing control.
            group->setCheckable(true);
            group->setChecked(toggle->isChecked());
            toggle->hide();
            QObject::connect(group, &QGroupBox::toggled, toggle, &QCheckBox::setChecked);
            QObject::connect(toggle, &QCheckBox::toggled, group, &QGroupBox::setChecked);
        }
        return group;
    }

    inline QGridLayout *grid(QGroupBox *group)
    {
        return static_cast<QGridLayout *>(group->layout());
    }

    inline void field(QGridLayout *layout, int row, const QString &title, QWidget *input, const QString &description)
    {
        auto *label = new QLabel(title);
        label->setBuddy(input);
        label->setToolTip(description);
        input->setToolTip(description);
        input->setAccessibleName(title);
        layout->addWidget(label, row, 0);
        layout->addWidget(input, row, 1);
        layout->setColumnStretch(1, 1);
    }
} // namespace model::controls
#endif
