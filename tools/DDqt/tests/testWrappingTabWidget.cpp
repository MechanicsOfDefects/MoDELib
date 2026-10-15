#include <WrappingTabWidget.h>
#include <DDqtControls.h>
#include <QLineEdit>
#include <QApplication>
#include <set>
#include <stdexcept>
#include <iostream>

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    try
    {
        QWidget panel;
        QCheckBox visibility(&panel);
        auto *group = model::controls::box("Test", &panel, "Group help", &visibility);
        auto *edit = new QLineEdit();
        model::controls::field(model::controls::grid(group), 0, "Value", edit, "Detailed field help");
        if(group->isChecked() || edit->isEnabled() || !visibility.isHidden())
        {
            throw std::runtime_error("Initial group state is incorrect");
        }
        group->setChecked(true);
        if(!visibility.isChecked() || !edit->isEnabled())
        {
            throw std::runtime_error("Group checkbox did not enable its controls");
        }
        visibility.setChecked(false);
        if(group->isChecked() || edit->isEnabled())
        {
            throw std::runtime_error("Programmatic visibility did not synchronize the group");
        }
        const auto *label = group->findChild<QLabel *>();
        if(!label || label->buddy() != edit || edit->toolTip().isEmpty() || label->toolTip().isEmpty())
        {
            throw std::runtime_error("Missing field label or tooltip");
        }
        auto *panelLayout = new QVBoxLayout(&panel);
        panelLayout->addWidget(group);
        auto *nested = new model::controls::GroupBox("Nested", group);
        nested->setCheckable(true);
        nested->setChecked(false);
        auto *nestedLayout = new QVBoxLayout(nested);
        auto *nestedEdit = new QLineEdit("kept value");
        nestedLayout->addWidget(nestedEdit);
        model::controls::grid(group)->addWidget(nested, 1, 0, 1, 2);
        auto *intentionallyHidden = new QLabel("Hidden option", group);
        model::controls::grid(group)->addWidget(intentionallyHidden, 2, 0);
        intentionallyHidden->hide();
        panel.show();
        app.processEvents();
        const int collapsedHeight = group->sizeHint().height();
        if(edit->isVisible())
        {
            throw std::runtime_error("Unchecked group contents are visible");
        }
        group->setChecked(true);
        app.processEvents();
        if(!edit->isVisible() || group->sizeHint().height() <= collapsedHeight || nestedEdit->isVisible())
        {
            throw std::runtime_error("Expansion failed or expanded an unchecked nested group");
        }
        nested->setChecked(true);
        group->setChecked(false);
        app.processEvents();
        if(edit->isVisible() || nestedEdit->isVisible() || group->sizeHint().height() != collapsedHeight)
        {
            throw std::runtime_error("Collapse did not reclaim space");
        }
        group->setChecked(true);
        app.processEvents();
        if(!nestedEdit->isVisible() || nestedEdit->text() != "kept value" || !intentionallyHidden->isHidden())
        {
            throw std::runtime_error("Collapse lost nested state, values or intentional visibility");
        }
        model::WrappingTabWidget tabs;
        for(const auto *name : {"Nodes", "Segments", "Loops", "Inclusions", "GlidePlanes", "Quadrature", "Chart",
                                "Fields", "SurfaceMesh", "SlipSystems"})
        {
            tabs.addTab(new QWidget(), name);
        }
        tabs.resize(340, 500);
        tabs.show();
        app.processEvents();
        auto buttons = tabs.findChildren<QPushButton *>();
        const auto checkLayout = [&]()
        {
            std::set<int> rows;
            for(const auto *button : buttons)
            {
                rows.insert(button->y());
                if(button->width() < button->sizeHint().width() ||
                   !button->parentWidget()->rect().contains(button->geometry()))
                {
                    throw std::runtime_error("Actor tab label is clipped");
                }
                for(const auto *other : buttons)
                {
                    if(button != other && button->geometry().intersects(other->geometry()))
                    {
                        throw std::runtime_error("Actor tab buttons overlap");
                    }
                }
            }
            return rows.size();
        };
        const auto narrowRows = checkLayout();
        if(narrowRows < 2)
        {
            throw std::runtime_error("Actor tabs did not wrap");
        }
        buttons[8]->click();
        if(tabs.currentIndex() != 8 || !buttons[8]->isChecked() || buttons[0]->isChecked())
        {
            throw std::runtime_error("Actor tab selection failed");
        }
        tabs.resize(1000, 500);
        app.processEvents();
        if(checkLayout() >= narrowRows || tabs.currentIndex() != 8)
        {
            throw std::runtime_error("Resize lost selection or failed to reflow");
        }
        tabs.resize(240, 500);
        app.processEvents();
        checkLayout();
        std::cout << "PASS: actor tab wrapping, full labels, selection and resize\n";
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
